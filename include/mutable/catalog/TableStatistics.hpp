#pragma once

#include <unordered_map>
#include <vector>
#include <string>
#include <memory>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <set>


namespace m
{

    namespace cnf {
        struct CNF;
    }
    namespace ast {
        struct Constant;
    }
    struct Table; // Forward declaration

    struct ColumnHistogram
    {
        // Equi-width histogram for numeric columns ONLY
        std::vector<std::size_t> bins;
        double min = 0.0;
        double max = 0.0;
        std::size_t num_distinct = 0;
        std::size_t null_count = 0;
        std::size_t total_count = 0;

        static const std::size_t DEFAULT_NUM_BINS = 10;

        // Default constructor
        ColumnHistogram() = default;

        // Create histogram from numeric values (only method needed)
        static ColumnHistogram create_numeric_histogram(const std::vector<double> &values,
                                                        std::size_t distinct_count,
                                                        std::size_t null_count,
                                                        std::size_t num_bins = DEFAULT_NUM_BINS);

        // Multiply two histograms for join estimation (assumes independence)
        ColumnHistogram multiply(const ColumnHistogram &other) const;

        // Merge two histograms (for union operations)
        ColumnHistogram merge(const ColumnHistogram &other) const;

        // Get selectivity for a range [low, high]
        double get_range_selectivity(double low, double high) const;

        // Get overall selectivity (distinct/total ratio)
        double get_selectivity() const;

        // Check if histogram is valid
        bool is_valid() const;

        ColumnHistogram filter_greater_than(double threshold) const;
        ColumnHistogram filter_less_than(double threshold) const;
        ColumnHistogram filter_range(double low, double high) const;
        ColumnHistogram filter_equal(double value) const;


    private:
        // Helper to align two histograms to same range and bin count
        static std::pair<ColumnHistogram, ColumnHistogram> align_histograms(
            const ColumnHistogram &left, const ColumnHistogram &right);
    };

    struct TableStatistics
    {
    public:
        std::string table_name;
        std::set<std::string> column_names;
        // Selectivity for ALL columns (table_name.column_name format)
        std::unordered_map<std::string, double> selectivity;
        // Histograms ONLY for numeric columns (table_name.column_name format)
        std::unordered_map<std::string, ColumnHistogram> histograms;
        std::size_t row_count = 0;

        // Compute statistics for a table
        void compute(const Table &table);

        // Extract table name from table
        void extract_table_name(const Table &table);

        // Get table name
        std::string get_table_name() const;

        void extract_column_names(const Table &table);

        // Get selectivity for any column (numeric or non-numeric)
        double get_selectivity(const std::string &table_col) const
        {
            auto it = selectivity.find(table_col);
            return it != selectivity.end() ? it->second : -1.0;
        }

        // Get histogram for numeric columns only (returns nullptr for non-numeric)
        const ColumnHistogram *get_histogram(const std::string &table_col) const
        {
            auto it = histograms.find(table_col);
            return it != histograms.end() ? &it->second : nullptr;
        }

        // Check if a column has a histogram (i.e., is numeric)
        bool has_histogram(const std::string &table_col) const
        {
            return histograms.find(table_col) != histograms.end();
        }

        // Multiply histograms for join estimation (only works for numeric columns)
        ColumnHistogram multiply_histograms(const std::string &left_col, const std::string &right_col) const;

        // Merge with another TableStatistics (for joins)
        TableStatistics merge_for_join(const TableStatistics &other) const;

        /**
         * Filter histograms based on CNF conditions
         * @param cnf_condition The CNF filter to apply
         * @return New TableStatistics with filtered histograms
         */
        TableStatistics filter_by_cnf(const cnf::CNF& cnf_condition) const;

    };
}; // namespace m