#pragma once

#include <unordered_map>
#include <vector>
#include <string>
#include <memory>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <mutable/util/Pool.hpp>

namespace m
{

    struct Table; // Forward declaration

    struct ColumnHistogram
    {
        // Simple equi-width histogram for numeric columns
        std::vector<std::size_t> bins;
        double min = 0.0;
        double max = 0.0;
        std::size_t num_distinct = 0;
        std::size_t null_count = 0;
    };

    struct TableStatistics
    {

        public:
        // Per-column selectivity (fraction of unique values)
        std::unordered_map<ThreadSafePooledString, double> selectivity;
        // Per-column histograms
        std::unordered_map<ThreadSafePooledString, ColumnHistogram> histograms;
        // Total row count
        std::size_t row_count = 0;

        // Compute statistics for a table (call after loading data)
        void compute(const Table &table);

        // Get selectivity for a column (fraction unique)
        double get_selectivity(const ThreadSafePooledString &col) const
        {
            auto it = selectivity.find(col);
            return it != selectivity.end() ? it->second : 1.0;
        }

        // Get histogram for a column
        const ColumnHistogram *get_histogram(const ThreadSafePooledString &col) const
        {
            auto it = histograms.find(col);
            return it != histograms.end() ? &it->second : nullptr;
        }

        // Overloads for std::string parameters (for convenience)
        double get_selectivity(const std::string &col) const
        {
            // Don't try to construct a ThreadSafePooledString directly
            // Instead, look through the keys for a match
            for (const auto &[key, val] : selectivity)
            {
                if (std::string(*key) == col)
                { // Convert the ThreadSafePooledString to std::string for comparison
                    return val;
                }
            }
            return -1.0;
        }

        const ColumnHistogram *get_histogram(const std::string &col) const
        {
            // Same approach for histograms
            for (const auto &[key, val] : histograms)
            {
                if (std::string(*key) == col)
                {
                    return &val;
                }
            }
            return nullptr;
        }
    };

} // namespace m