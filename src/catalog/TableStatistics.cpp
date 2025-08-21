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

        // FIXED: Cartesian product within each bin
        for (std::size_t i = 0; i < result.bins.size(); ++i)
        {
            // For equi-join: multiply bin counts directly (cartesian product within bin)
            if (aligned_left.bins[i] > 0 && aligned_right.bins[i] > 0)
            {
                // Check for potential overflow before multiplication
                if (aligned_left.bins[i] > std::numeric_limits<std::size_t>::max() / aligned_right.bins[i])
                {
                    // Use conservative estimate to prevent overflow
                    result.bins[i] = std::min(aligned_left.bins[i], aligned_right.bins[i]);
                }
                else
                {
                    result.bins[i] = aligned_left.bins[i] * aligned_right.bins[i];
                }
            }
            else
            {
                result.bins[i] = 0;
            }
        }

        // Conservative join selectivity: min of distinct counts
        result.num_distinct = std::min(aligned_left.num_distinct, aligned_right.num_distinct);
        result.null_count = 0; // Nulls don't participate in joins
        result.total_count = std::accumulate(result.bins.begin(), result.bins.end(), std::size_t(0));

        // SAFETY CHECK: Cap result at reasonable bounds
        std::size_t max_possible_join =
            (aligned_left.total_count - aligned_left.null_count) *
            (aligned_right.total_count - aligned_right.null_count);
        if (result.total_count > max_possible_join)
        {
            // Scale down proportionally if result is too large
            double scale_factor = double(max_possible_join) / double(result.total_count);
            for (auto &bin : result.bins)
            {
                bin = static_cast<std::size_t>(bin * scale_factor);
            }
            result.total_count = std::accumulate(result.bins.begin(), result.bins.end(), std::size_t(0));
        }

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
        if (bins.empty())
            return false;
        if (total_count == 0)
            return false;
        if (min > max)
            return false;
        if (num_distinct > total_count)
            return false;
        if (null_count > total_count)
            return false;

        // Check that bin counts sum correctly
        std::size_t bin_sum = std::accumulate(bins.begin(), bins.end(), std::size_t(0));
        if (bin_sum + null_count != total_count)
            return false;

        // Check for overflow in bins
        for (auto bin_count : bins)
        {
            if (bin_count > total_count)
                return false;
        }

        return true;
    }

    std::pair<ColumnHistogram, ColumnHistogram> ColumnHistogram::align_histograms(
        const ColumnHistogram &left, const ColumnHistogram &right)
    {
        if (left.bins.empty() || right.bins.empty())
        {
            return {left, right}; // Return as-is for empty histograms
        }

        // Find common range and bin count
        double common_min = std::min(left.min, right.min);
        double common_max = std::max(left.max, right.max);
        std::size_t common_bins = std::max(left.bins.size(), right.bins.size());

        // If histograms are already aligned, return as-is
        if (left.min == right.min && left.max == right.max &&
            left.bins.size() == right.bins.size())
        {
            return {left, right};
        }

        // Create aligned versions with uniform redistribution
        ColumnHistogram aligned_left = redistribute_to_range(left, common_min, common_max, common_bins);
        ColumnHistogram aligned_right = redistribute_to_range(right, common_min, common_max, common_bins);

        return {aligned_left, aligned_right};
    }

    // Helper method for uniform redistribution
    ColumnHistogram ColumnHistogram::redistribute_to_range(
        const ColumnHistogram &hist, double new_min, double new_max, std::size_t new_bin_count)
    {
        ColumnHistogram result;
        result.min = new_min;
        result.max = new_max;
        result.bins.resize(new_bin_count, 0);
        result.null_count = hist.null_count;
        result.num_distinct = hist.num_distinct;
        result.total_count = hist.total_count;

        // If original histogram is empty or degenerate, return empty result
        if (hist.bins.empty() || hist.max <= hist.min || new_max <= new_min)
        {
            return result;
        }

        // Calculate bin widths
        double old_bin_width = (hist.max - hist.min) / hist.bins.size();
        double new_bin_width = (new_max - new_min) / new_bin_count;

        // Redistribute each old bin uniformly across new bins
        for (std::size_t old_bin = 0; old_bin < hist.bins.size(); ++old_bin)
        {
            if (hist.bins[old_bin] == 0)
                continue;

            double old_bin_start = hist.min + old_bin * old_bin_width;
            double old_bin_end = hist.min + (old_bin + 1) * old_bin_width;

            // Find overlap with new range
            double overlap_start = std::max(old_bin_start, new_min);
            double overlap_end = std::min(old_bin_end, new_max);

            if (overlap_start >= overlap_end)
                continue; // No overlap

            // Distribute this old bin's count across overlapping new bins
            for (std::size_t new_bin = 0; new_bin < new_bin_count; ++new_bin)
            {
                double new_bin_start = new_min + new_bin * new_bin_width;
                double new_bin_end = new_min + (new_bin + 1) * new_bin_width;

                // Calculate overlap between old bin and new bin
                double bin_overlap_start = std::max(overlap_start, new_bin_start);
                double bin_overlap_end = std::min(overlap_end, new_bin_end);

                if (bin_overlap_start < bin_overlap_end)
                {
                    // Calculate fraction of old bin that goes into this new bin
                    double overlap_width = bin_overlap_end - bin_overlap_start;
                    double old_bin_kept_width = overlap_end - overlap_start;
                    double fraction = overlap_width / old_bin_kept_width;

                    result.bins[new_bin] += static_cast<std::size_t>(hist.bins[old_bin] * fraction);
                }
            }
        }

        return result;
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

        std::vector<std::vector<std::string>> string_values(schema.num_entries());
        std::vector<std::vector<double>> numeric_values(schema.num_entries());
        std::vector<bool> is_numeric(schema.num_entries(), false);
        std::vector<std::size_t> null_counts(schema.num_entries(), 0);

        // Prepare NDV/selectivity computation using Tuple
        for (std::size_t col = 0; col < schema.num_entries(); ++col)
        {
            std::string col_name = std::string(*schema[col].id.name);
            std::string full_key = table_name + "." + col_name;

            std::unordered_map<Value, int> value_count;

            // For each row, extract the value as a single-column Tuple
            auto loader = Interpreter::compile_load(schema, table.store().memory().addr(), table.layout(), schema, 0, 0);
            Tuple tuple(schema);

            for (std::size_t row = 0; row < row_count; ++row)
            {
                
                Tuple *args[] = {&tuple};
                loader(args);

                if (!tuple.is_null(col))
                {
                    ++value_count[tuple[col]];
                }
            }

            std::size_t nd = value_count.size();
            distinct_counts[full_key] = nd;

            // Compute Top K TODO: Make this a setting
            
            int most_frequent_value_count = -1;

            for (const auto& [value, count] : value_count){
                if (count > most_frequent_value_count){
                    most_frequent_value_count = count;
                }
            }
            
            most_frequent_values[full_key] = most_frequent_value_count;

            double sel = row_count > 0 ? double(nd) / double(row_count) : 1.0;
            selectivity[full_key] = sel;

            if (is_numeric[col] && !numeric_values[col].empty())
            {
                histograms[full_key] = ColumnHistogram::create_numeric_histogram(
                    numeric_values[col], nd, null_counts[col]);
            }
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
        // combine names, columns, histograms, selectivity as before…
        result.table_name = table_name + "_" + other.table_name;
        result.column_names = column_names;
        result.column_names.insert(other.column_names.begin(), other.column_names.end());

        result.selectivity = selectivity;
        for (auto &kv : other.selectivity)
            result.selectivity[kv.first] = kv.second;

        result.histograms = histograms;
        for (auto &kv : other.histograms)
            result.histograms[kv.first] = kv.second;

        // merge distinct_counts too
        result.distinct_counts = distinct_counts;
        for (auto &kv : other.distinct_counts)
            result.distinct_counts[kv.first] = kv.second;

        // row_count will be set by the estimator
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

        // Step 1: Create new histogram with adjusted range [threshold, max]
        ColumnHistogram result;
        result.min = threshold;             // Adjust min to threshold
        result.max = max;                   // Keep original max
        result.bins.resize(bins.size(), 0); // Same number of bins, all zero
        result.null_count = null_count;

        // Calculate old and new bin widths
        double old_bin_width = (max - min) / bins.size();
        double new_bin_width = (result.max - result.min) / result.bins.size();

        // Step 2: For each NEW bin, calculate overlap with OLD bins
        for (std::size_t new_bin = 0; new_bin < result.bins.size(); ++new_bin)
        {
            double new_bin_start = result.min + new_bin * new_bin_width;
            double new_bin_end = result.min + (new_bin + 1) * new_bin_width;

            std::size_t new_bin_count = 0;

            // Check overlap with each OLD bin
            for (std::size_t old_bin = 0; old_bin < bins.size(); ++old_bin)
            {
                double old_bin_start = min + old_bin * old_bin_width;
                double old_bin_end = min + (old_bin + 1) * old_bin_width;

                // Only consider the part of old bin that's > threshold
                double effective_old_start = std::max(old_bin_start, threshold);
                if (old_bin_end <= threshold)
                {
                    continue; // This old bin is completely below threshold
                }

                // Calculate overlap between new bin and (filtered) old bin
                double overlap_start = std::max(new_bin_start, effective_old_start);
                double overlap_end = std::min(new_bin_end, old_bin_end);

                if (overlap_start < overlap_end)
                {
                    // There's overlap - calculate what fraction of the old bin's count to add
                    double old_bin_kept_width = old_bin_end - effective_old_start;
                    double overlap_width = overlap_end - overlap_start;
                    double fraction = overlap_width / old_bin_kept_width;

                    // Calculate how much of the old bin was kept (due to filtering)
                    double old_bin_total_width = old_bin_end - old_bin_start;
                    if (old_bin_total_width <= 0)
                    {
                        continue; // Skip degenerate bins
                    }
                    double kept_fraction = old_bin_kept_width / old_bin_total_width;
                    std::size_t old_bin_kept_count = static_cast<std::size_t>(bins[old_bin] * kept_fraction);

                    // Add the overlapping portion
                    new_bin_count += static_cast<std::size_t>(old_bin_kept_count * fraction);
                }
            }

            result.bins[new_bin] = new_bin_count;
        }

        // Step 3: Update metadata
        result.total_count = std::accumulate(result.bins.begin(), result.bins.end(), std::size_t(0)) + result.null_count;

        // Scale distinct count by fraction of data kept
        double fraction = double(result.total_count - result.null_count) / double(total_count - null_count);
        result.num_distinct = static_cast<std::size_t>(num_distinct * fraction);

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

        // Step 1: Create new histogram with adjusted range [min, threshold]
        ColumnHistogram result;
        result.min = min;                   // Keep original min
        result.max = threshold;             // Adjust max to threshold
        result.bins.resize(bins.size(), 0); // Same number of bins, all zero
        result.null_count = null_count;

        // Calculate old and new bin widths
        double old_bin_width = (max - min) / bins.size();
        double new_bin_width = (result.max - result.min) / result.bins.size();

        // Step 2: For each NEW bin, calculate overlap with OLD bins
        for (std::size_t new_bin = 0; new_bin < result.bins.size(); ++new_bin)
        {
            double new_bin_start = result.min + new_bin * new_bin_width;
            double new_bin_end = result.min + (new_bin + 1) * new_bin_width;

            std::size_t new_bin_count = 0;

            // Check overlap with each OLD bin
            for (std::size_t old_bin = 0; old_bin < bins.size(); ++old_bin)
            {
                double old_bin_start = min + old_bin * old_bin_width;
                double old_bin_end = min + (old_bin + 1) * old_bin_width;

                // Only consider the part of old bin that's < threshold
                double effective_old_end = std::min(old_bin_end, threshold);
                if (old_bin_start >= threshold)
                {
                    continue; // This old bin is completely above threshold
                }

                // Calculate overlap between new bin and (filtered) old bin
                double overlap_start = std::max(new_bin_start, old_bin_start);
                double overlap_end = std::min(new_bin_end, effective_old_end);

                if (overlap_start < overlap_end)
                {
                    // There's overlap - calculate what fraction of the old bin's count to add
                    double old_bin_kept_width = effective_old_end - old_bin_start;
                    double overlap_width = overlap_end - overlap_start;
                    double fraction = overlap_width / old_bin_kept_width;

                    // Calculate how much of the old bin was kept (due to filtering)
                    double old_bin_total_width = old_bin_end - old_bin_start;
                    if (old_bin_total_width <= 0)
                    {
                        continue; // Skip degenerate bins
                    }
                    double kept_fraction = old_bin_kept_width / old_bin_total_width;
                    std::size_t old_bin_kept_count = static_cast<std::size_t>(bins[old_bin] * kept_fraction);

                    // Add the overlapping portion
                    new_bin_count += static_cast<std::size_t>(old_bin_kept_count * fraction);
                }
            }

            result.bins[new_bin] = new_bin_count;
        }

        // Step 3: Update metadata
        result.total_count = std::accumulate(result.bins.begin(), result.bins.end(), std::size_t(0)) + result.null_count;

        // Scale distinct count by fraction of data kept
        double fraction = double(result.total_count - result.null_count) / double(total_count - null_count);
        result.num_distinct = static_cast<std::size_t>(num_distinct * fraction);

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

        // Step 1: Create new histogram with single bin for the exact value
        ColumnHistogram result;
        result.min = value;       // Min equals the target value
        result.max = value;       // Max equals the target value
        result.bins.resize(1, 0); // Single bin, start with 0
        result.null_count = null_count;
        result.num_distinct = 1; // Only one distinct value

        // Step 2: Find which old bin contains this value and apply uniformity assumption
        double old_bin_width = (max - min) / bins.size();
        std::size_t target_bin = static_cast<std::size_t>((value - min) / old_bin_width);
        target_bin = std::min(target_bin, bins.size() - 1); // Clamp to valid range

        // Apply uniformity assumption: assume values are evenly distributed within the bin
        if (bins[target_bin] > 0 && num_distinct > 0)
        {
            // Calculate how many distinct values are approximately in this bin
            double distinct_per_bin = double(num_distinct) / double(bins.size());

            // If this bin likely contains multiple distinct values,
            // assume uniform distribution and take 1/distinct_per_bin fraction
            if (distinct_per_bin > 1.0)
            {
                double fraction = 1.0 / distinct_per_bin;
                result.bins[0] = static_cast<std::size_t>(bins[target_bin] * fraction);
            }
            else
            {
                // If distinct_per_bin <= 1, this bin might contain only one distinct value
                // So we can take all of it (conservative estimate)
                result.bins[0] = bins[target_bin];
            }
        }
        else
        {
            result.bins[0] = 0; // No data in the target bin
        }

        // Step 3: Update metadata
        result.total_count = result.bins[0] + result.null_count;

        return result;
    }

    /**
     * Rescales all histograms to match a specified target cardinality.
     *
     * @param target_cardinality The desired cardinality to scale all histograms to
     * @return A new TableStatistics object with rescaled histograms
     */
    TableStatistics TableStatistics::rescale_histograms_to_cardinality(std::size_t target_cardinality) const
    {
        // Create a copy of the current statistics
        TableStatistics result = *this;

        // If there are no histograms or target is 0, nothing to do
        if (histograms.empty() || target_cardinality == 0)
        {
            return result;
        }

        // Scale each histogram directly to the target cardinality
        for (auto &[col_name, hist] : result.histograms)
        {
            if (!hist.is_valid() || hist.total_count == 0 || hist.total_count == target_cardinality)
            {
                continue; // Skip invalid or already-at-target histograms
            }

            // Calculate scaling factor directly for this histogram
            double hist_scale_factor = double(target_cardinality) / double(hist.total_count);

            // Scale non-null counts proportionally
            std::size_t non_null = hist.total_count - hist.null_count;
            std::size_t new_non_null = static_cast<std::size_t>(double(non_null) * hist_scale_factor);

            // Scale each bin proportionally
            if (non_null > 0)
            {
                double bin_scale = double(new_non_null) / double(non_null);
                for (auto &bin : hist.bins)
                {
                    bin = static_cast<std::size_t>(double(bin) * bin_scale);
                }
            }

            // Update total count while preserving null count
            hist.total_count = new_non_null + hist.null_count;
        }

        // Update the table row count to match the target cardinality
        result.row_count = target_cardinality;

        // Update distinct counts differently depending on scaling direction
        bool scaling_up = target_cardinality > row_count;
        for (auto &[col_name, distinct_count] : result.distinct_counts)
        {
            if (scaling_up)
            {
                // When scaling up, conservatively increase NDV
                // Scale by sqrt of the ratio to model diminishing new distinct values
                double scale_ratio = sqrt(double(target_cardinality) / double(row_count));
                distinct_count = std::min(target_cardinality,
                                          static_cast<std::size_t>(distinct_count * scale_ratio));
            }
            else if (distinct_count > target_cardinality)
            {
                // When scaling down, cap at new cardinality
                distinct_count = target_cardinality;
            }
        }

        // Update selectivity values based on new cardinality and distinct counts
        for (auto &[col_name, sel] : result.selectivity)
        {
            if (result.distinct_counts.count(col_name) && result.row_count > 0)
            {
                sel = double(result.distinct_counts.at(col_name)) / double(result.row_count);
            }
        }

        return result;
    }

    // TODO: This function also has to rescale every other histogram according to the card count
    TableStatistics TableStatistics::filter_by_cnf(const cnf::CNF &cnf_condition) const
    {
        TableStatistics result = *this;

        for (const auto &clause : cnf_condition)
        {
            if (clause.size() != 1)
            {
                continue;
            }

            const auto &predicate = clause[0];

            if (predicate.negative())
            {
                continue;
            }

            if (auto binary_expr = cast<const ast::BinaryExpr>(&predicate.expr()))
            {

                auto lhs = cast<const ast::Designator>(binary_expr->lhs.get());
                if (!lhs || !lhs->has_table_name())
                {
                    continue;
                }

                std::string table_col = std::string(*lhs->table_name.text) + "." +
                                        std::string(*lhs->attr_name.text);

                // Check if we have a histogram for this column
                auto hist_it = result.histograms.find(table_col);
                if (hist_it == result.histograms.end())
                {
                    continue;
                }

                // Extract right side (should be a constant value) - SIMPLE WAY
                if (auto constant = cast<const ast::Constant>(binary_expr->rhs.get()))
                {

                    std::ostringstream oss;
                    oss << *constant;
                    std::string value_str = oss.str();

                    double filter_value;
                    try
                    {
                        filter_value = std::stod(value_str); 
                    }
                    catch (const std::exception &)
                    {
                        continue;
                    }

                    ColumnHistogram filtered_hist;
                    switch (binary_expr->op().type)
                    {
                    case TK_LESS:
                        filtered_hist = hist_it->second.filter_less_than(filter_value);
                        break;
                    case TK_LESS_EQUAL:
                        filtered_hist = hist_it->second.filter_less_than(filter_value);
                        break;
                    case TK_GREATER:
                        filtered_hist = hist_it->second.filter_greater_than(filter_value);
                        break;
                    case TK_GREATER_EQUAL:
                        filtered_hist = hist_it->second.filter_greater_than(filter_value);
                        break;
                    case TK_EQUAL:
                        filtered_hist = hist_it->second.filter_equal(filter_value);
                        break;
                    default:
                        continue; // Unsupported operator
                    }

                    hist_it->second = filtered_hist;
                }
            }
        }

        std::size_t min_cardinality = std::numeric_limits<std::size_t>::max();
        for (const auto &[col_name, hist] : result.histograms)
        {
            if (hist.is_valid() && hist.total_count > 0 && hist.total_count < min_cardinality)
            {
                min_cardinality = hist.total_count;
            }
        }

        // Only rescale if we found valid histograms
        if (min_cardinality != std::numeric_limits<std::size_t>::max())
        {
            return result.rescale_histograms_to_cardinality(min_cardinality);
        }

        // If no valid histograms found, return as-is
        return result;
    }

    std::size_t TableStatistics::estimate_group_by_cardinality(const std::vector<std::string> &group_columns) const
    {
        if (group_columns.empty())
        {
            return row_count; // No grouping
        }

        // For multiple columns, multiply distinct counts (assuming independence)
        std::size_t estimated_groups = 1;
        for (const std::string &col : group_columns)
        {
            const auto *hist = get_histogram(col);
            if (hist && hist->is_valid())
            {
                estimated_groups *= hist->num_distinct;
            }
            else
            {
                // Fallback: assume high cardinality for non-numeric columns
                estimated_groups *= row_count;
            }
        }

        return estimated_groups;
    }

    TableStatistics TableStatistics::apply_group_by(const std::vector<std::string> &group_columns) const
    {
        TableStatistics result = *this;

        // Update row count to estimated group count
        result.row_count = estimate_group_by_cardinality(group_columns);

        // Apply group by to all histograms
        for (auto &[col_name, histogram] : result.histograms)
        {
            if (std::find(group_columns.begin(), group_columns.end(), col_name) != group_columns.end())
            {
                // This column is being grouped by
                histogram = histogram.apply_group_by();
            }
            else
            {
                // This column is not being grouped by - estimate its new distribution
                // For simplicity, assume it maintains the same distinct count but fewer total values
                double scale_factor = double(result.row_count) / double(row_count);
                histogram.total_count = static_cast<std::size_t>(histogram.total_count * scale_factor);
                for (auto &bin : histogram.bins)
                {
                    bin = static_cast<std::size_t>(bin * scale_factor);
                }
            }
        }

        return result;
    }

} // namespace m