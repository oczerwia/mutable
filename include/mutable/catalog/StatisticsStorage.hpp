#pragma once

#include <unordered_map>
#include <string>
#include <memory>
#include <iostream>
#include <mutable/catalog/TableStatistics.hpp>

namespace m
{

    class StatisticsStorage
    {
    public:
        // Singleton access
        static StatisticsStorage &Get()
        {
            static StatisticsStorage instance;
            return instance;
        }

        // Store or update statistics for a table
        void set_statistics(const std::string &table_name, std::unique_ptr<TableStatistics> stats)
        {
            table_stats_[table_name] = std::move(stats);
        }

        // Get statistics for a table (nullptr if not found)
        TableStatistics *get_statistics(const std::string &table_name)
        {
            auto it = table_stats_.find(table_name);
            return it != table_stats_.end() ? it->second.get() : nullptr;
        }

        // Print all selectivities for all tables
        void print_all_selectivities() const
        {
            for (const auto &[table, stats_ptr] : table_stats_)
            {
                if (stats_ptr)
                {
                    for (const auto &[col, sel] : stats_ptr->selectivity)
                    {
                        std::cout << "Selectivity for table '" << table
                                  << "', column '" << *col << "': " << sel << "\n";
                    }
                }
            }
        }

    private:
        StatisticsStorage() = default;
        std::unordered_map<std::string, std::unique_ptr<TableStatistics>> table_stats_;
    };

} // namespace m