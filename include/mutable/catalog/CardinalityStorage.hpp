#pragma once

#include <mutable/IR/PlanTable.hpp>
#include <mutable/IR/Operator.hpp>
#include <mutable/IR/CNF.hpp>
#include <mutable/util/ADT.hpp>

#include <unordered_map>
#include <memory>
#include <functional>
#include <utility>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <set>

namespace m
{
    /**
     * @brief Complete data structure containing all relevant cardinality information for easy extraction
     */
    struct CardinalityData
    {
        // Cardinality information
        double estimated_cardinality = -1.0; // Estimated by optimizer
        double true_cardinality = -1.0;      // Actual from execution

        std::pair<double,double> estimated_range = {-1.0, -1.0};

        // Subproblem identification
        Subproblem subproblem; // Complete bitmask

        // Tables involved
        std::vector<std::string> table_names; // Table names in this subproblem

        // Operation information
        bool has_filter = false; // Has filter operation
        bool has_join = false;   // Has join operation

        // Performance metrics
        double selectivity = -1.0;   // Computed selectivity if available
        double error_percent = -1.0; // Estimation error percentage

        // NEW: Store filters as strings for this specific subproblem
        std::set<std::string> filter_strings;
    };

    struct ReducedQueryGraph

    // add a function with which you can compare the reducedquerygraph with a set of table names or build the reduced query graph from the query graph in the beginning
    // then iterate over the list of sets and return the cardinality if it is there.
    {
    public:
        // needed data

        // CNF::filter;

        std::set<std::string> source_names_;
        std::set<std::string> source_filters_;
        double cardinality_;

        // NEW: Add range field
        std::pair<double, double> cardinality_range_ = {-1.0, -1.0}; // -1 indicates unset

        bool has_range() const
        {
            return cardinality_range_.first >= 0.0; // Range exists if set to valid values
        }

        void set_range(const std::pair<double, double> &range)
        {
            cardinality_range_ = range;
            // Keep point estimate (backward compatibility)
            // cardinality_ already set separately
        }

        std::pair<double, double> get_range() const
        {
            return cardinality_range_;
        }

        // comparison data
        // sources already defined above

        // filter
        // joins
        // group_by
        // aggregate

        bool operator==(const ReducedQueryGraph &other) const
        {
            return source_names_ == other.source_names_ &&
                   source_filters_ == other.source_filters_;
        }

        bool operator!=(const ReducedQueryGraph &other) const
        {
            return !(*this == other);
        }
        double get_cardinality() const
        {
            return this->cardinality_;
        }
    };

    /**
     * @brief Singleton class that stores cardinality information from query execution
     */
    class CardinalityStorage
    {
    private:
        // Efficiently store data with shared_ptr to avoid duplication
        std::vector<std::shared_ptr<CardinalityData>> all_cardinality_data_;

        // Map from subproblem to data index for fast lookup
        std::unordered_map<Subproblem, std::size_t, SubproblemHash> subproblem_to_data_;

        bool debug_output_ = true;

        // NEW TRY
        std::vector<std::string> current_table_names_;
        std::set<std::string> current_filters_;
        // List to store reduced query graphs for cardinality lookup
        std::vector<ReducedQueryGraph> reduced_query_graphs_;
        // when we find a cardinality, we store it on this variable
        double cardinality_;

        // Storage for found range (similar to cardinality_)
        std::pair<double, double> cardinality_range_ = {-1.0, -1.0};

    public:
        // Delete copy/move constructors and assignment operators
        CardinalityStorage(const CardinalityStorage &) = delete;
        CardinalityStorage &operator=(const CardinalityStorage &) = delete;
        CardinalityStorage(CardinalityStorage &&) = delete;
        CardinalityStorage &operator=(CardinalityStorage &&) = delete;

        /**
         * @brief Get the singleton instance
         */
        static CardinalityStorage &Get()
        {
            static CardinalityStorage instance;
            return instance;
        }

        void quick_test_generator()
        {
            ReducedQueryGraph test_query;
            test_query.source_names_.insert("table_1");
            test_query.source_names_.insert("table_2");
            test_query.source_names_.insert("table_3");
            test_query.cardinality_ = 10.0;

            reduced_query_graphs_.push_back(test_query);

            ReducedQueryGraph test_query_2;
            test_query_2.source_names_.insert("table_1");
            test_query_2.source_names_.insert("table_2");
            test_query_2.cardinality_ = 20.0;

            reduced_query_graphs_.push_back(test_query_2);

            // Add a third test query with same tables as test_query_2 but with a filter
            ReducedQueryGraph test_query_3;
            test_query_3.source_names_.insert("table_1");
            test_query_3.source_names_.insert("table_2");
            test_query_3.source_filters_.insert("(table_1.col_1 < 400)");
            test_query_3.cardinality_ = 15.0;

            reduced_query_graphs_.push_back(test_query_3);

            if (debug_output_)
            {
                std::cout << "Generated test query with tables: ";
                for (const auto &table : test_query.source_names_)
                {
                    std::cout << table << " ";
                }
                std::cout << "and cardinality: " << test_query.cardinality_ << std::endl;
            }
        }

        /**
         * @brief Get the current table names
         *
         * @return const std::set<std::string>& Reference to current table names
         */
        const std::vector<std::string> &get_current_table_names() const
        {
            return current_table_names_;
        }

        /**
         * @brief Update the current table names set by traversing the plan table
         *
         * @param plan_table The plan table to extract table names from
         */
        void update_current_table_names(const QueryGraph &query_graph)
        {
            std::vector<std::string> table_names;

            // Traverse all data sources in the query graph
            const auto &sources = query_graph.sources();
            for (const auto &source_ptr : sources)
            {
                const DataSource &ds = *source_ptr;
                auto name_opt = ds.name();
                if (name_opt.has_value())
                {
                    table_names.push_back(*name_opt);
                }
            }

            // Update the current table names vector
            current_table_names_ = table_names;

            if (debug_output_)
            {
                std::cout << "Updated current table names from QueryGraph: ";
                for (const auto &name : current_table_names_)
                {
                    std::cout << name << " ";
                }
                std::cout << std::endl;
            }
        }

        bool has_stored_cardinality(SmallBitset involved_tables)
        {
            // 1. Iterate over the Bitset
            // 2. get the involved tables from previously determined set of involved tables
            // 3. iterate over the vector of all previously stored queries and search for equal
            //      sets of tables

            std::set<std::string> table_names;
            for (auto pos : involved_tables)
            {
                if (pos < current_table_names_.size())
                {
                    table_names.insert(current_table_names_[pos]);
                }
            }
            for (const auto &stored_query : reduced_query_graphs_)
            {
                if (stored_query.source_names_ == table_names)
                {
                    if (stored_query.source_filters_ == current_filters_)
                    {
                        cardinality_ = stored_query.get_cardinality();
                        return true;
                    }
                }
            }
            return false;
        }

        /**
         * @brief Get debug output status
         */
        bool debug_output() const
        {
            return debug_output_;
        }

        /**
         * @brief Get the stored cardinality value
         */
        double get_cardinality() const
        {
            return cardinality_;
        }

        void extract_all_filters_as_strings(const QueryGraph &G)
        {
            std::set<std::string> filter_strings;
            const auto &sources = G.sources();
            for (const auto &source_ptr : sources)
            {
                const DataSource &ds = *source_ptr;
                if (!ds.filter().empty())
                {
                    // Convert to string and normalize (sort clauses/predicates)
                    std::string filter_str = to_string(ds.filter());
                    filter_strings.insert(filter_str);
                }
            }
            current_filters_ = filter_strings;
        }

        /**
         * @brief Extract metadata from operators (tables, filters, joins) - SIMPLIFIED
         *
         * @param op The operator to extract metadata from
         * @param data The CardinalityData to update
         */
        void extract_operator_metadata_(const Operator &op, CardinalityData &data)
        {
            if (debug_output_)
            {
                std::cout << "  extract_operator_metadata_: subproblem = " << data.subproblem << std::endl;
                std::cout << "  current_table_names_ size = " << current_table_names_.size() << std::endl;
            }

            // Extract table names directly from the subproblem bitmask using current_table_names_
            for (unsigned pos = 0; pos < current_table_names_.size() && pos < 64; ++pos)
            {
                if (data.subproblem[pos])
                {
                    const std::string &table_name = current_table_names_[pos];

                    if (debug_output_)
                    {
                        std::cout << "    Found table at pos " << pos << ": " << table_name << std::endl;
                    }

                    // Only add if not already present
                    if (std::find(data.table_names.begin(), data.table_names.end(), table_name) == data.table_names.end())
                    {
                        data.table_names.push_back(table_name);
                    }
                }
            }

            // ONLY extract direct filters from FilterOperator (no complex propagation)
            if (auto filter_op = dynamic_cast<const FilterOperator *>(&op))
            {
                data.has_filter = true;
                std::ostringstream ss;
                ss << filter_op->filter();
                std::string filter_str = ss.str();
                data.filter_strings.insert(filter_str);

                if (debug_output_)
                {
                    std::cout << "  Stored direct filter: " << filter_str << std::endl;
                }
            }

            // Mark join operations
            if (auto join_op = dynamic_cast<const JoinOperator *>(&op))
            {
                data.has_join = true;
            }
        }

        /**
         * @brief Traverses the operator tree and collects cardinality information
         *
         * @param op The current operator being processed
         */
        void traverse_operator_tree_(const Operator &op)
        {
            if (debug_output_)
            {
                std::cout << "Processing operator: " << typeid(op).name() << std::endl;
            }

            if (op.has_info())
            {
                const Subproblem &subproblem = op.info().subproblem;

                // Create or find the cardinality data for this subproblem
                std::shared_ptr<CardinalityData> data;
                auto it = subproblem_to_data_.find(subproblem);

                if (it == subproblem_to_data_.end())
                {
                    // Create new data entry
                    data = std::make_shared<CardinalityData>();
                    data->subproblem = subproblem;

                    // Add to our collections
                    subproblem_to_data_[subproblem] = all_cardinality_data_.size();
                    all_cardinality_data_.push_back(data);
                }
                else
                {
                    // Use existing entry
                    data = all_cardinality_data_[it->second];
                }

                // Update cardinality information
                if (op.has_info() && op.info().estimated_cardinality > 0)
                {
                    data->estimated_cardinality = op.info().estimated_cardinality;

                    if (op.info().estimated_range.first >= 0.0)
                    {
                        data->estimated_range = op.info().estimated_range;
                    }
                }
                data->true_cardinality = op.get_emitted_tuples();

                // Calculate error percentage if both values are available
                if (data->estimated_cardinality > 0 && data->true_cardinality > 0)
                {
                    data->error_percent = std::abs(data->estimated_cardinality - data->true_cardinality) /
                                          data->true_cardinality * 100.0;
                }

                // Calculate selectivity if processed count is available
                size_t processed = op.get_processed_tuples();
                if (processed > 0)
                {
                    data->selectivity = static_cast<double>(data->true_cardinality) / processed;
                }

                // Extract metadata from operator
                extract_operator_metadata_(op, *data);
            }

            // Recursively traverse children
            if (auto consumer = dynamic_cast<const Consumer *>(&op))
            {
                for (auto child : consumer->children())
                {
                    traverse_operator_tree_(*child);
                }
            }
        }

        void map_true_cardinalities_to_logical_plan_(const Operator &root)
        {
            if (debug_output_)
            {
                std::cout << "Starting to map true cardinalities to logical plan..." << std::endl;
            }

            // Clear previous data
            all_cardinality_data_.clear();
            subproblem_to_data_.clear();

            // Traverse and collect cardinalities
            traverse_operator_tree_(root);

            // Debug output
            if (debug_output_)
            {
                std::cout << "\nCollected " << all_cardinality_data_.size() << " cardinality entries." << std::endl;
                for (const auto &data : all_cardinality_data_)
                {
                    std::cout << "  Subproblem: estimated=" << data->estimated_cardinality;
                    if (data->estimated_range.first >= 0.0)
                    {
                        std::cout << ", range_estimated=[" << data->estimated_range.first 
                                << "-" << data->estimated_range.second << "]";
                    }
                    
                    std::cout << ", actual=" << data->true_cardinality << std::endl;

                    std::cout << "    Tables: ";
                    for (const auto &name : data->table_names)
                    {
                        std::cout << name << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }

        /**
         * @brief Extract true cardinalities and populate reduced query graphs
         *
         * @param root The root operator of the executed physical plan
         */
        void extract_reduced_query_graphs_from_execution(const Operator &root)
        {
            if (debug_output_)
            {
                std::cout << "Extracting reduced query graphs from execution..." << std::endl;
            }

            // First, run the new cardinality extraction
            map_true_cardinalities_to_logical_plan_(root);

            // Step 1: Create temporary vector of all new reduced query graphs
            std::vector<ReducedQueryGraph> temp_graphs;

            for (const auto &cardinality_data : all_cardinality_data_)
            {
                if (cardinality_data->true_cardinality > 0)
                {
                    ReducedQueryGraph new_graph;

                    // Set tables
                    for (const auto &table_name : cardinality_data->table_names)
                    {
                        new_graph.source_names_.insert(table_name);
                    }

                    // Set ONLY direct filters from this cardinality_data
                    for (const auto &filter_str : cardinality_data->filter_strings)
                    {
                        new_graph.source_filters_.insert(filter_str);
                    }

                    new_graph.cardinality_ = cardinality_data->true_cardinality;
                    if (cardinality_data->estimated_range.first >= 0.0)
                        new_graph.set_range(cardinality_data->estimated_range);
                    temp_graphs.push_back(new_graph);
                }
            }

            // Step 2: Propagate filters to multi-table entries
            for (auto &graph : temp_graphs)
            {
                if (graph.source_names_.size() > 1) // Multi-table graphs
                {
                    // Check each table in this graph
                    for (const auto &table_name : graph.source_names_)
                    {
                        // Look through all single-table graphs for filters that apply to this table
                        for (const auto &single_graph : temp_graphs)
                        {
                            if (single_graph.source_names_.size() == 1 &&
                                single_graph.source_names_.count(table_name) > 0 &&
                                !single_graph.source_filters_.empty())
                            {
                                // Add all filters from this single-table graph
                                for (const auto &filter_str : single_graph.source_filters_)
                                {
                                    graph.source_filters_.insert(filter_str);

                                    if (debug_output_)
                                    {
                                        std::cout << "  Propagated filter '" << filter_str
                                                  << "' to multi-table graph containing " << table_name << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Step 3: Update or add to storage
            for (const auto &new_graph : temp_graphs)
            {
                bool found_existing = false;
                for (auto &existing_graph : reduced_query_graphs_)
                {
                    if (existing_graph.source_names_ == new_graph.source_names_ &&
                        existing_graph.source_filters_ == new_graph.source_filters_)
                    {
                        existing_graph.cardinality_ = new_graph.cardinality_;
                        if (existing_graph.cardinality_range_.first >= 0.0){
                            existing_graph.cardinality_range_ = new_graph.cardinality_range_;
                        }
                        found_existing = true;
                        if (debug_output_)
                        {
                            std::cout << "  Updated existing ReducedQueryGraph: tables={";
                            for (const auto &table : existing_graph.source_names_)
                            {
                                std::cout << table << " ";
                            }
                            std::cout << "}, filters={";
                            for (const auto &filter : existing_graph.source_filters_)
                            {
                                std::cout << filter << " ";
                            }
                            std::cout << "}, cardinality=" << existing_graph.cardinality_;
                            // Add range output
                            if (existing_graph.cardinality_range_.first >= 0.0)
                            {
                                std::cout << ", range=[" << existing_graph.cardinality_range_.first
                                        << "-" << existing_graph.cardinality_range_.second << "]";
                            }
                            std::cout << std::endl;
                        }
                        break;
                    }
                }

                // Only add if not found existing
                if (!found_existing)
                {
                    reduced_query_graphs_.push_back(new_graph);

                    if (debug_output_)
                    {
                        std::cout << "  Added new ReducedQueryGraph: tables={";
                        for (const auto &table : new_graph.source_names_)
                        {
                            std::cout << table << " ";
                        }
                        std::cout << "}, filters={";
                        for (const auto &filter : new_graph.source_filters_)
                        {
                            std::cout << filter << " ";
                        }
                        std::cout << "}, cardinality=" << new_graph.cardinality_ << std::endl;
                    }
                }
            }

            if (debug_output_)
            {
                std::cout << "Total ReducedQueryGraphs in storage: " << reduced_query_graphs_.size() << std::endl;
            }
        }

        // NEW: Check if a range exists
        bool has_stored_cardinality_range(SmallBitset involved_tables)
        {
            // Reuse existing logic to find the correct ReducedQueryGraph
            if (has_stored_cardinality(involved_tables))
            {
                // Find the graph we just matched in has_stored_cardinality
                std::set<std::string> table_names;
                for (auto pos : involved_tables)
                {
                    if (pos < current_table_names_.size())
                    {
                        table_names.insert(current_table_names_[pos]);
                    }
                }

                for (const auto &stored_query : reduced_query_graphs_)
                {
                    if (stored_query.source_names_ == table_names &&
                        stored_query.source_filters_ == current_filters_)
                    {
                        if (stored_query.has_range())
                        {
                            cardinality_range_ = stored_query.get_range();
                            return true;
                        }
                        return false;
                    }
                }
            }
            return false;
        }

        std::pair<double, double> get_cardinality_range() const
        {
            return cardinality_range_;
        }

        void store_cardinality_range(const std::set<std::string> &tables,
                                     const std::pair<double, double> &range,
                                     double true_cardinality)
        {
            // Find or create the ReducedQueryGraph
            for (auto &graph : reduced_query_graphs_)
            {
                if (graph.source_names_ == tables &&
                    graph.source_filters_ == current_filters_)
                {
                    // Update existing graph
                    graph.cardinality_ = true_cardinality;
                    graph.set_range(range);
                    return;
                }
            }

            // Create new graph if not found
            ReducedQueryGraph new_graph;
            new_graph.source_names_ = tables;
            new_graph.source_filters_ = current_filters_;
            new_graph.cardinality_ = true_cardinality;
            new_graph.set_range(range);
            reduced_query_graphs_.push_back(new_graph);
        }

    private:
        // Add default constructor
        CardinalityStorage() = default;
    };
} // namespace m
