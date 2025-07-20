#include <mutable/catalog/TableStatistics.hpp>
#include <mutable/catalog/Schema.hpp>
#include <mutable/IR/Tuple.hpp>
#include <backend/Interpreter.hpp>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <limits>
#include <numeric>

namespace m
{
    // ColumnHistogram implementations (simplified for numeric-only)
    ColumnHistogram ColumnHistogram::create_numeric_histogram(const std::vector<double> &values,
                                                              std::size_t distinct_count,
                                                              std::size_t null_count,
                                                              std::size_t num_bins)
    {
        ColumnHistogram hist;
        hist.num_distinct = distinct_count;
        hist.null_count = null_count;
        hist.total_count = values.size() + null_count;

        if (values.empty())
        {
            return hist;
        }

        // Find min and max
        auto minmax = std::minmax_element(values.begin(), values.end());
        hist.min = *minmax.first;
        hist.max = *minmax.second;

        // Initialize bins
        hist.bins.resize(num_bins, 0);

        // Create equi-width bins
        if (hist.max > hist.min)
        {
            double bin_width = (hist.max - hist.min) / num_bins;

            for (double value : values)
            {
                std::size_t bin_index = static_cast<std::size_t>((value - hist.min) / bin_width);
                if (bin_index >= num_bins)
                    bin_index = num_bins - 1;
                hist.bins[bin_index]++;
            }
        }
        else
        {
            // All values are the same
            hist.bins[0] = values.size();
        }

        return hist;
    }

    ColumnHistogram ColumnHistogram::multiply(const ColumnHistogram &other) const
    {
        if (bins.empty() || other.bins.empty())
        {
            return ColumnHistogram(); // Return empty histogram
        }

        // Align histograms to same range and bin count
        auto [aligned_left, aligned_right] = align_histograms(*this, other);

        ColumnHistogram result;
        result.min = std::max(aligned_left.min, aligned_right.min);
        result.max = std::min(aligned_left.max, aligned_right.max);
        result.bins.resize(aligned_left.bins.size());

        // Multiply corresponding bins (assuming independence)
        double left_total = aligned_left.total_count - aligned_left.null_count;
        double right_total = aligned_right.total_count - aligned_right.null_count;

        for (std::size_t i = 0; i < result.bins.size(); ++i)
        {
            double left_prob = left_total > 0 ? double(aligned_left.bins[i]) / left_total : 0.0;
            double right_prob = right_total > 0 ? double(aligned_right.bins[i]) / right_total : 0.0;
            result.bins[i] = static_cast<std::size_t>(left_prob * right_prob * left_total * right_total);
        }

        result.num_distinct = std::min(aligned_left.num_distinct, aligned_right.num_distinct);
        result.null_count = 0; // Nulls don't participate in joins
        result.total_count = std::accumulate(result.bins.begin(), result.bins.end(), std::size_t(0));

        return result;
    }

    ColumnHistogram ColumnHistogram::merge(const ColumnHistogram &other) const
    {
        if (bins.empty() || other.bins.empty())
        {
            return ColumnHistogram(); // Return empty histogram
        }

        // Align histograms
        auto [aligned_left, aligned_right] = align_histograms(*this, other);

        ColumnHistogram result;
        result.min = std::min(aligned_left.min, aligned_right.min);
        result.max = std::max(aligned_left.max, aligned_right.max);
        result.bins.resize(aligned_left.bins.size());

        // Add corresponding bins
        for (std::size_t i = 0; i < result.bins.size(); ++i)
        {
            result.bins[i] = aligned_left.bins[i] + aligned_right.bins[i];
        }

        result.num_distinct = num_distinct + other.num_distinct;
        result.null_count = null_count + other.null_count;
        result.total_count = total_count + other.total_count;

        return result;
    }

    double ColumnHistogram::get_range_selectivity(double low, double high) const
    {
        if (bins.empty() || total_count == 0)
        {
            return 0.1; // Default selectivity
        }

        if (max <= min)
        {
            return (low <= min && min <= high) ? 1.0 : 0.0;
        }

        double bin_width = (max - min) / bins.size();
        std::size_t start_bin = static_cast<std::size_t>((low - min) / bin_width);
        std::size_t end_bin = static_cast<std::size_t>((high - min) / bin_width);

        start_bin = std::min(start_bin, bins.size() - 1);
        end_bin = std::min(end_bin, bins.size() - 1);

        std::size_t count = 0;
        for (std::size_t i = start_bin; i <= end_bin; ++i)
        {
            count += bins[i];
        }

        return double(count) / double(total_count - null_count);
    }

    double ColumnHistogram::get_selectivity() const
    {
        if (total_count == 0)
            return 0.0;
        return double(num_distinct) / double(total_count);
    }

    bool ColumnHistogram::is_valid() const
    {
        return total_count > 0 && !bins.empty();
    }

    std::pair<ColumnHistogram, ColumnHistogram> ColumnHistogram::align_histograms(
        const ColumnHistogram &left, const ColumnHistogram &right)
    {

        if (left.bins.empty() || right.bins.empty())
        {
            return {left, right}; // Return as-is for empty histograms
        }

        // Find common range
        double common_min = std::min(left.min, right.min);
        double common_max = std::max(left.max, right.max);
        std::size_t common_bins = std::max(left.bins.size(), right.bins.size());

        // Create aligned histograms
        ColumnHistogram aligned_left = left;
        ColumnHistogram aligned_right = right;

        // Ensure same bin count and range
        aligned_left.min = common_min;
        aligned_left.max = common_max;
        aligned_right.min = common_min;
        aligned_right.max = common_max;

        if (aligned_left.bins.size() != common_bins)
        {
            aligned_left.bins.resize(common_bins, 0);
        }
        if (aligned_right.bins.size() != common_bins)
        {
            aligned_right.bins.resize(common_bins, 0);
        }

        return {aligned_left, aligned_right};
    }

    // TableStatistics implementations
    void TableStatistics::compute(const Table &table)
    {
        row_count = table.store().num_rows();
        selectivity.clear();
        histograms.clear();

        extract_table_name(table);
        extract_column_names(table);

        const auto &schema = table.schema();
        Tuple tuple(schema);

        // Prepare containers for each column
        std::vector<std::vector<std::string>> string_values(schema.num_entries());
        std::vector<std::vector<double>> numeric_values(schema.num_entries());
        std::vector<bool> is_numeric(schema.num_entries(), false);
        std::vector<std::size_t> null_counts(schema.num_entries(), 0);

        // For each row, collect values
        for (std::size_t row = 0; row < row_count; ++row)
        {
            auto loader = Interpreter::compile_load(schema, table.store().memory().addr(), table.layout(), schema, row, 0);
            Tuple *args[] = {&tuple};
            loader(args);

            for (std::size_t col = 0; col < schema.num_entries(); ++col)
            {
                auto &entry = tuple[col];

                if (tuple.is_null(col))
                {
                    null_counts[col]++;
                    continue;
                }

                std::ostringstream oss;
                oss << entry;
                std::string str_value = oss.str();
                string_values[col].push_back(str_value);

                // Check if numeric and collect numeric values
                const auto &type = schema[col].type;
                if (type->is_integral() || type->is_floating_point())
                {
                    try
                    {
                        double numeric_val = std::stod(str_value);
                        numeric_values[col].push_back(numeric_val);
                        is_numeric[col] = true;
                    }
                    catch (const std::exception &)
                    {
                        // Conversion failed, treat as non-numeric
                        is_numeric[col] = false;
                    }
                }
            }
        }

        // Compute statistics for each column
        for (std::size_t col = 0; col < schema.num_entries(); ++col)
        {
            std::string col_name = std::string(*schema[col].id.name);
            std::string full_key = table_name + "." + col_name;

            // Compute selectivity for ALL columns (numeric and non-numeric)
            std::unordered_set<std::string> unique_vals(string_values[col].begin(), string_values[col].end());
            double sel = row_count > 0 ? double(unique_vals.size()) / double(row_count) : 1.0;
            selectivity[full_key] = sel;

            // Create histogram ONLY for numeric columns
            if (is_numeric[col] && !numeric_values[col].empty())
            {
                histograms[full_key] = ColumnHistogram::create_numeric_histogram(
                    numeric_values[col], unique_vals.size(), null_counts[col]);
            }
            // Note: No histogram created for non-numeric columns
        }
    }

    ColumnHistogram TableStatistics::multiply_histograms(const std::string &left_col, const std::string &right_col) const
    {
        auto left_hist = get_histogram(left_col);
        auto right_hist = get_histogram(right_col);

        if (!left_hist || !right_hist)
        {
            return ColumnHistogram(); // Return empty histogram if either column is not numeric
        }

        return left_hist->multiply(*right_hist);
    }

    TableStatistics TableStatistics::merge_for_join(const TableStatistics &other) const
    {
        TableStatistics result;

        // Combine table names
        result.table_name = table_name + "_" + other.table_name;

        // Merge column names
        result.column_names = column_names;
        result.column_names.insert(other.column_names.begin(), other.column_names.end());

        // Merge selectivities (for all columns)
        result.selectivity = selectivity;
        for (const auto &[key, value] : other.selectivity)
        {
            result.selectivity[key] = value;
        }

        // Merge histograms (only numeric columns will have histograms)
        result.histograms = histograms;
        for (const auto &[key, value] : other.histograms)
        {
            result.histograms[key] = value;
        }

        // Row count will be set by join estimation
        result.row_count = 0;

        return result;
    }

    std::string TableStatistics::get_table_name() const
    {
        return table_name;
    }

    void TableStatistics::extract_table_name(const Table &table)
    {
        table_name = std::string(*table.name());
    }

    void TableStatistics::extract_column_names(const Table &table)
    {
        column_names.clear();
        const auto &schema = table.schema();

        for (std::size_t col = 0; col < schema.num_entries(); ++col)
        {
            std::string col_name = std::string(*schema[col].id.name);
            column_names.insert(col_name);
        }
    }

    ColumnHistogram ColumnHistogram::filter_greater_than(double threshold) const
    {
        if (bins.empty() || total_count == 0 || threshold >= max)
        {
            return ColumnHistogram(); // Empty result
        }

        if (threshold <= min)
        {
            return *this; // No filtering needed
        }

        ColumnHistogram result = *this;
        double bin_width = (max - min) / bins.size();

        // Find the first bin that contains values > threshold
        std::size_t start_bin = static_cast<std::size_t>((threshold - min) / bin_width);
        start_bin = std::min(start_bin, bins.size() - 1);

        // Zero out bins before threshold
        for (std::size_t i = 0; i < start_bin; ++i)
        {
            result.bins[i] = 0;
        }

        // Partially zero out the bin containing the threshold
        if (start_bin < bins.size())
        {
            double bin_start = min + start_bin * bin_width;
            double bin_end = min + (start_bin + 1) * bin_width;

            if (threshold > bin_start && threshold < bin_end)
            {
                // Assume uniform distribution within the bin
                double fraction_kept = (bin_end - threshold) / (bin_end - bin_start);
                result.bins[start_bin] = static_cast<std::size_t>(result.bins[start_bin] * fraction_kept);
            }
        }

        // Update metadata
        result.min = threshold;
        result.total_count = std::accumulate(result.bins.begin(), result.bins.end(), std::size_t(0)) + result.null_count;

        // Estimate new distinct count (conservative)
        double selectivity_ratio = double(result.total_count - result.null_count) / double(total_count - null_count);
        result.num_distinct = static_cast<std::size_t>(num_distinct * selectivity_ratio);

        return result;
    }

    ColumnHistogram ColumnHistogram::filter_less_than(double threshold) const
    {
        if (bins.empty() || total_count == 0 || threshold <= min)
        {
            return ColumnHistogram(); // Empty result
        }

        if (threshold >= max)
        {
            return *this; // No filtering needed
        }

        ColumnHistogram result = *this;
        double bin_width = (max - min) / bins.size();

        // Find the last bin that contains values < threshold
        std::size_t end_bin = static_cast<std::size_t>((threshold - min) / bin_width);
        end_bin = std::min(end_bin, bins.size() - 1);

        // Zero out bins after threshold
        for (std::size_t i = end_bin + 1; i < bins.size(); ++i)
        {
            result.bins[i] = 0;
        }

        // Partially zero out the bin containing the threshold
        if (end_bin < bins.size())
        {
            double bin_start = min + end_bin * bin_width;
            double bin_end = min + (end_bin + 1) * bin_width;

            if (threshold > bin_start && threshold < bin_end)
            {
                // Assume uniform distribution within the bin
                double fraction_kept = (threshold - bin_start) / (bin_end - bin_start);
                result.bins[end_bin] = static_cast<std::size_t>(result.bins[end_bin] * fraction_kept);
            }
        }

        // Update metadata
        result.max = threshold;
        result.total_count = std::accumulate(result.bins.begin(), result.bins.end(), std::size_t(0)) + result.null_count;

        // Estimate new distinct count
        double selectivity_ratio = double(result.total_count - result.null_count) / double(total_count - null_count);
        result.num_distinct = static_cast<std::size_t>(num_distinct * selectivity_ratio);

        return result;
    }

    ColumnHistogram ColumnHistogram::filter_range(double low, double high) const
    {
        return filter_greater_than(low).filter_less_than(high);
    }

    ColumnHistogram ColumnHistogram::filter_equal(double value) const
    {
        if (bins.empty() || total_count == 0 || value < min || value > max)
        {
            return ColumnHistogram(); // Empty result
        }

        ColumnHistogram result;
        result.min = value;
        result.max = value;
        result.bins.resize(1);
        result.null_count = null_count;

        // Find which bin contains this value
        double bin_width = (max - min) / bins.size();
        std::size_t target_bin = static_cast<std::size_t>((value - min) / bin_width);
        target_bin = std::min(target_bin, bins.size() - 1);

        // Estimate count for this specific value
        // Assume uniform distribution within the bin
        if (num_distinct > 0)
        {
            double avg_count_per_distinct = double(total_count - null_count) / double(num_distinct);
            result.bins[0] = static_cast<std::size_t>(avg_count_per_distinct);
        }
        else
        {
            result.bins[0] = 0;
        }

        result.total_count = result.bins[0] + result.null_count;
        result.num_distinct = result.total_count > result.null_count ? 1 : 0;

        return result;
    }
} // namespace m