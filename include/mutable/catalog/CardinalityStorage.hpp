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
    // Simple filter predicate representation that doesn't rely on specific expression types
    struct FilterPredicate
    {
        std::string table_name;
        std::string column_name;
        std::string op;    // "=", ">", "<", ">=", "<=", "<>", etc.
        std::string value; // String representation of the value

        // Equality comparison for exact matching
        bool operator==(const FilterPredicate &other) const
        {
            return table_name == other.table_name &&
                   column_name == other.column_name &&
                   op == other.op &&
                   value == other.value;
        }

        std::string to_string() const
        {
            return table_name + "." + column_name + " " + op + " " + value;
        }
    };

    /**
     * @brief Complete data structure containing all relevant cardinality information for easy extraction
     */
    struct CardinalityData
    {
        // Cardinality information
        double estimated_cardinality = -1.0; // Estimated by optimizer
        double true_cardinality = -1.0;      // Actual from execution

        // Subproblem identification
        Subproblem subproblem; // Complete bitmask

        // Tables involved
        std::vector<std::string> table_names; // Table names in this subproblem

        // Table position mapping (which bit position corresponds to which table)
        std::vector<std::pair<unsigned, std::string>> table_positions; // Bit position -> table name

        // Operation information
        bool has_filter = false; // Has filter operation
        bool has_join = false;   // Has join operation

        // Performance metrics
        double selectivity = -1.0;   // Computed selectivity if available
        double error_percent = -1.0; // Estimation error percentage

        // Add filter information
        std::vector<FilterPredicate> filter_predicates;

        // NEW: Store filters as strings for this specific subproblem
        std::set<std::string> filter_strings;

        // Get the full bitmask as a string (e.g., "1101")
        std::string get_bitmask_string(size_t max_positions) const
        {
            std::string result(max_positions, '0');
            for (auto &[pos, _] : table_positions)
            {
                if (pos < max_positions)
                    result[pos] = '1';
            }
            return result;
        }

        // Get a table name by its bit position
        std::string get_table_name_by_position(unsigned position) const
        {
            for (const auto &[pos, name] : table_positions)
            {
                if (pos == position)
                    return name;
            }
            return "";
        }

        // Check if a specific table is included
        bool includes_table(const std::string &table_name) const
        {
            for (const auto &name : table_names)
            {
                if (name == table_name)
                    return true;
            }
            return false;
        }

        // Helper to check if filter predicates match
        bool filter_predicates_match(const std::vector<FilterPredicate> &other_filters) const
        {
            if (filter_predicates.size() != other_filters.size())
                return false;

            // Check each predicate (order independent)
            for (const auto &pred : filter_predicates)
            {
                bool found = false;
                for (const auto &other : other_filters)
                {
                    if (pred == other)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    return false;
            }

            return true;
        }
    };

    /**
     * @brief Custom data object to attach to PlanTableEntry
     */
    struct PlanTableEntryCardinalityData : public PlanTableEntryData
    {
        std::shared_ptr<CardinalityData> data; // Use shared_ptr for efficient copying

        PlanTableEntryCardinalityData(std::shared_ptr<CardinalityData> data_ptr) : data(data_ptr) {}
    };

    struct StoredQueryPlan
    {
        // Simplified representation of the plan structure
        std::unordered_map<Subproblem, Subproblem, SubproblemHash> join_structure; // maps subproblem -> left child

        // Add filter predicates
        std::vector<FilterPredicate> filter_predicates;

        // Map from subproblem to the filters that apply to it
        std::unordered_map<Subproblem, std::vector<size_t>, SubproblemHash> subproblem_filters; // Maps to indices in filter_predicates

        // Map from subproblem to the tables it contains
        std::unordered_map<Subproblem, std::vector<std::string>, SubproblemHash> subproblem_tables;
    };

    class ReducedQueryGraph

    // add a function with which you can compare the reducedquerygraph with a set of table names or build the reduced query graph from the query graph in the beginning
    // then iterate over the list of sets and return the cardinality if it is there.
    {
    public:
        // needed data

        // CNF::filter;

        std::set<std::string> source_names_;
        std::set<std::string> source_filters_;
        double cardinality_;

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

        // Global table position mapping
        std::vector<std::pair<unsigned, std::string>> global_table_positions_;
        size_t max_bit_positions_ = 0;

        bool debug_output_ = true;

        // Add to the private section of CardinalityStorage class
        std::vector<StoredQueryPlan> stored_query_plans_;

        // Private constructor for singleton pattern
        CardinalityStorage() = default;

        // NEW TRY
        std::vector<std::string> current_table_names_;
        std::set<std::string> current_filters_;
        // List to store reduced query graphs for cardinality lookup
        std::vector<ReducedQueryGraph> reduced_query_graphs_;
        // when we find a cardinality, we store it on this variable
        double cardinality_;

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

        double get_cardinality() const
        {
            return this->cardinality_;
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

            // Note: We no longer need collect_table_positions since we use current_table_names_

            // Traverse and collect cardinalities
            traverse_operator_tree_(root);

            // Debug output
            if (debug_output_)
            {
                std::cout << "\nCollected " << all_cardinality_data_.size() << " cardinality entries." << std::endl;
                for (const auto &data : all_cardinality_data_)
                {
                    std::cout << "  Subproblem: estimated=" << data->estimated_cardinality
                              << ", actual=" << data->true_cardinality << std::endl;

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
                // Check if this combination already exists
                bool found_existing = false;
                for (auto &existing_graph : reduced_query_graphs_)
                {
                    if (existing_graph.source_names_ == new_graph.source_names_ &&
                        existing_graph.source_filters_ == new_graph.source_filters_)
                    {
                        // Update existing entry
                        existing_graph.cardinality_ = new_graph.cardinality_;
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
                            std::cout << "}, cardinality=" << existing_graph.cardinality_ << std::endl;
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

        /**
         * @brief Maps actual cardinalities from physical operators back to logical plan subproblems
         *
         * @param root The root operator of the executed physical plan
         */
        void map_true_cardinalities_to_logical_plan(const Operator &root)
        {
            std::cout << "Starting to map true cardinalities to logical plan..." << std::endl;

            // Clear previous data
            all_cardinality_data_.clear();
            subproblem_to_data_.clear();
            global_table_positions_.clear();

            // First pass: collect all table names and assign bit positions
            collect_table_positions(root);

            // Second pass: traverse and collect cardinalities
            traverse_operator_tree(root);

            // Debug output
            if (debug_output_)
            {
                std::cout << "\nTable to bit position mapping:" << std::endl;
                for (const auto &[pos, table] : global_table_positions_)
                {
                    std::cout << "  Bit " << pos << " = " << table << std::endl;
                }

                std::cout << "\nCollected " << all_cardinality_data_.size() << " cardinality entries." << std::endl;
                for (const auto &data : all_cardinality_data_)
                {
                    std::cout << "  Subproblem [" << data->get_bitmask_string(max_bit_positions_)
                              << "]: estimated=" << data->estimated_cardinality
                              << ", actual=" << data->true_cardinality << std::endl;

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
         * @brief First pass: collect all table names and assign bit positions
         *
         * @param root The root operator
         */
        void collect_table_positions(const Operator &root)
        {
            unsigned next_position = 0;

            std::function<void(const Operator &)> traverse = [&](const Operator &op)
            {
                if (auto scan_op = dynamic_cast<const ScanOperator *>(&op))
                {
                    try
                    {
                        std::ostringstream ss;
                        ss << scan_op->alias();
                        std::string table_name = ss.str();

                        if (!table_name.empty())
                        {
                            // Check if we already have a position for this table
                            bool found = false;
                            for (const auto &[pos, name] : global_table_positions_)
                            {
                                if (name == table_name)
                                {
                                    found = true;
                                    break;
                                }
                            }

                            // If not found, assign the next available position
                            if (!found)
                            {
                                global_table_positions_.emplace_back(next_position, table_name);
                                next_position++;
                                max_bit_positions_ = next_position; // Update max positions

                                if (debug_output_)
                                {
                                    std::cout << "  Assigned bit position " << (next_position - 1)
                                              << " to table " << table_name << std::endl;
                                }
                            }
                        }
                    }
                    catch (...)
                    {
                        // Ignore failures
                    }
                }

                // Recursively traverse children
                if (auto consumer = dynamic_cast<const Consumer *>(&op))
                {
                    for (auto child : consumer->children())
                    {
                        traverse(*child);
                    }
                }
            };

            traverse(root);

            if (debug_output_)
            {
                std::cout << "Collected " << global_table_positions_.size() << " tables." << std::endl;
            }
        }

        /**
         * @brief Traverses the operator tree and collects cardinality information
         *
         * @param op The current operator being processed
         */
        void traverse_operator_tree(const Operator &op)
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
                extract_operator_metadata(op, *data);

                // Update table positions based on the subproblem bitmask
                update_table_positions(*data);
            }

            // Recursively traverse children
            if (auto consumer = dynamic_cast<const Consumer *>(&op))
            {
                for (auto child : consumer->children())
                {
                    traverse_operator_tree(*child);
                }
            }
        }

        /**
         * @brief Update table positions for a CardinalityData based on its subproblem
         *
         * @param data The CardinalityData to update
         */
        void update_table_positions(CardinalityData &data)
        {
            // Clear existing positions
            data.table_positions.clear();

            // For each potential bit position (using the full range)
            for (unsigned pos = 0; pos < max_bit_positions_; ++pos)
            {
                // If this bit is set in the subproblem
                if (data.subproblem[pos])
                {
                    // Look up the table name for this position
                    for (const auto &[map_pos, table_name] : global_table_positions_)
                    {
                        if (map_pos == pos)
                        {
                            data.table_positions.emplace_back(pos, table_name);
                            break;
                        }
                    }
                }
            }
        }

        /**
         * @brief Extract metadata from operators (tables, filters, joins)
         *
         * @param op The operator to extract metadata from
         * @param data The CardinalityData to update
         */
        void extract_operator_metadata(const Operator &op, CardinalityData &data)
        {
            // For ScanOperator - extract table name
            if (auto scan_op = dynamic_cast<const ScanOperator *>(&op))
            {
                try
                {
                    std::ostringstream ss;
                    ss << scan_op->alias();
                    std::string table_name = ss.str();

                    if (!table_name.empty())
                    {
                        // Only add if not already present
                        if (std::find(data.table_names.begin(), data.table_names.end(), table_name) == data.table_names.end())
                        {
                            data.table_names.push_back(table_name);
                        }
                    }
                }
                catch (...)
                {
                    // Silent failure
                }
            }

            // Enhanced FilterOperator handling with string-based parsing
            if (auto filter_op = dynamic_cast<const FilterOperator *>(&op))
            {
                data.has_filter = true;

                // Get the filter CNF in string form
                std::ostringstream ss;
                ss << filter_op->filter();
                std::string filter_str = ss.str();

                if (debug_output_)
                {
                    std::cout << "  Extracting filter: " << filter_str << std::endl;
                }

                // First, split by AND if there are multiple clauses
                std::vector<std::string> clauses;
                size_t start = 0;
                size_t pos;
                while ((pos = filter_str.find(" AND ", start)) != std::string::npos)
                {
                    clauses.push_back(filter_str.substr(start, pos - start));
                    start = pos + 5; // 5 is the length of " AND "
                }
                clauses.push_back(filter_str.substr(start)); // Add the last clause

                // Process each clause
                for (const auto &clause : clauses)
                {
                    // Extract filter components using common operators - order matters for proper parsing
                    std::vector<std::string> operators = {">=", "<=", "<>", "!=", "=", ">", "<"};

                    // Find which operator is in this clause
                    size_t op_pos = std::string::npos;
                    std::string found_op;
                    for (const auto &op : operators)
                    {
                        size_t pos = clause.find(op);
                        if (pos != std::string::npos)
                        {
                            op_pos = pos;
                            found_op = op;
                            break;
                        }
                    }

                    if (op_pos != std::string::npos)
                    {
                        // Found an operator, extract left and right sides
                        std::string left = clause.substr(0, op_pos);
                        std::string right = clause.substr(op_pos + found_op.length());

                        // Trim whitespace
                        auto trim = [](std::string &s)
                        {
                            s.erase(0, s.find_first_not_of(" \t"));
                            s.erase(s.find_last_not_of(" \t") + 1);
                        };

                        trim(left);
                        trim(right);

                        // Parse left side to extract table and column
                        size_t dot_pos = left.find('.');
                        if (dot_pos != std::string::npos)
                        {
                            std::string table = left.substr(0, dot_pos);
                            std::string column = left.substr(dot_pos + 1);

                            // Create filter predicate
                            FilterPredicate pred;
                            pred.table_name = table;
                            pred.column_name = column;
                            pred.op = found_op;
                            pred.value = right;

                            if (debug_output_)
                            {
                                std::cout << "    Extracted predicate: " << pred.to_string() << std::endl;
                            }

                            data.filter_predicates.push_back(pred);
                        }
                        else if (debug_output_)
                        {
                            std::cout << "    Could not parse predicate (no table.column format): " << clause << std::endl;
                        }
                    }
                    else if (debug_output_)
                    {
                        std::cout << "    Could not find operator in: " << clause << std::endl;
                    }
                }

                if (debug_output_)
                {
                    std::cout << "  Extracted " << data.filter_predicates.size() << " predicates from filter" << std::endl;
                }
            }

            // For JoinOperator - mark presence and gather source tables
            if (auto join_op = dynamic_cast<const JoinOperator *>(&op))
            {
                data.has_join = true;

                // For joins, merge source tables from children
                if (auto consumer = dynamic_cast<const Consumer *>(&op))
                {
                    for (auto child : consumer->children())
                    {
                        if (child->has_info())
                        {
                            auto it = subproblem_to_data_.find(child->info().subproblem);
                            if (it != subproblem_to_data_.end())
                            {
                                const auto &child_data = all_cardinality_data_[it->second];

                                // Add source tables from child to this operator's info
                                for (const auto &table : child_data->table_names)
                                {
                                    if (std::find(data.table_names.begin(), data.table_names.end(), table) == data.table_names.end())
                                    {
                                        data.table_names.push_back(table);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /**
         * @brief Returns a logical tree representation with both estimated and true cardinalities
         *
         * @param PT The plan table representing the logical plan
         * @return std::unique_ptr<PlanTable> The augmented logical plan
         */
        template <typename PlanTable>
        std::unique_ptr<PlanTable> get_logical_tree_with_cardinalities(const PlanTable &original_PT)
        {
            std::cout << "Creating logical tree with cardinalities..." << std::endl;
            auto PT_copy = std::make_unique<PlanTable>(original_PT);

            // For each entry in the plan table
            for (std::size_t i = 1; i < PT_copy->size(); ++i)
            {
                Subproblem s(i);
                if (PT_copy->has_plan(s))
                {
                    auto &entry = (*PT_copy)[s];

                    // If we have cardinality info for this subproblem
                    auto it = subproblem_to_data_.find(s);
                    if (it != subproblem_to_data_.end())
                    {
                        // Attach the data to the plan table entry
                        entry.data = std::make_unique<PlanTableEntryCardinalityData>(
                            all_cardinality_data_[it->second]);

                        if (debug_output_)
                        {
                            const auto &data = all_cardinality_data_[it->second];
                            std::cout << "  Added cardinality data to subproblem ["
                                      << data->get_bitmask_string(max_bit_positions_) << "]" << std::endl;
                        }
                    }
                }
            }

            std::cout << "Logical tree with cardinalities created successfully." << std::endl;
            return PT_copy;
        }

        /**
         * @brief Get all collected cardinality data
         *
         * @return const std::vector<std::shared_ptr<CardinalityData>>& All cardinality data
         */
        const std::vector<std::shared_ptr<CardinalityData>> &get_all_cardinality_data() const
        {
            return all_cardinality_data_;
        }

        /**
         * @brief Get cardinality data for a specific subproblem
         *
         * @param subproblem The subproblem to lookup
         * @return std::shared_ptr<CardinalityData> The data (or nullptr if not found)
         */
        std::shared_ptr<CardinalityData> get_cardinality_data(const Subproblem &subproblem) const
        {
            auto it = subproblem_to_data_.find(subproblem);
            if (it != subproblem_to_data_.end())
            {
                return all_cardinality_data_[it->second];
            }
            return nullptr;
        }

        /**
         * @brief Get tables that appear in all cardinality data entries
         *
         * @return std::vector<std::string> List of all table names
         */
        std::vector<std::string> get_all_tables() const
        {
            std::vector<std::string> result;
            for (const auto &[pos, name] : global_table_positions_)
            {
                result.push_back(name);
            }
            return result;
        }

        /**
         * @brief Toggle debug output
         */
        void set_debug_output(bool enable)
        {
            debug_output_ = enable;
        }

        bool debug_output() const
        {
            return debug_output_;
        }

        /**
         * @brief Store a query plan for future matching
         *
         * @param plan_table The plan table containing the query plan
         */
        template <typename PlanTable>
        void store_query_plan(const PlanTable &plan_table)
        {
            StoredQueryPlan stored_plan;

            // Extract the join structure from the plan table
            for (std::size_t i = 1; i < plan_table.size(); ++i)
            {
                Subproblem s(i);
                if (plan_table.has_plan(s))
                {
                    const auto &entry = plan_table[s];

                    // Store join structure
                    if (s.size() > 1 && entry.left && entry.right)
                    {
                        stored_plan.join_structure[s] = entry.left;
                    }

                    // Extract and store filter predicates if present
                    if (entry.data)
                    {
                        if (auto cardinality_data = dynamic_cast<PlanTableEntryCardinalityData *>(entry.data.get()))
                        {
                            if (cardinality_data->data && cardinality_data->data->has_filter)
                            {
                                // Add filter predicates
                                for (const auto &filter_pred : cardinality_data->data->filter_predicates)
                                {
                                    // Store filter predicate
                                    size_t pred_idx = stored_plan.filter_predicates.size();
                                    stored_plan.filter_predicates.push_back(filter_pred);

                                    // Map the subproblem to this filter predicate
                                    stored_plan.subproblem_filters[s].push_back(pred_idx);
                                }
                            }

                            // Store table names for this subproblem
                            if (cardinality_data->data)
                            {
                                stored_plan.subproblem_tables[s] = cardinality_data->data->table_names;

                                if (debug_output_)
                                {
                                    std::cout << "  Stored tables for subproblem " << s << ": ";
                                    for (const auto &table : cardinality_data->data->table_names)
                                        std::cout << table << " ";
                                    std::cout << std::endl;
                                }
                            }
                        }
                    }
                }
            }

            stored_query_plans_.push_back(std::move(stored_plan));

            if (debug_output_)
            {
                std::cout << "Stored query plan #" << stored_query_plans_.size()
                          << " with " << stored_plan.join_structure.size() << " joins"
                          << " and " << stored_plan.filter_predicates.size() << " filters" << std::endl;
            }
        }

        /**
         * @brief Find stored query plans that match the given plan
         *
         * @param plan_table The plan table to match against stored plans
         * @param exact_match If true, requires exact join structure and filter match
         * @return std::vector<size_t> Indices of matching stored plans
         */
        template <typename PlanTable>
        std::vector<size_t> find_matching_query_plans(const PlanTable &plan_table, bool exact_match = true) const
        {
            std::vector<size_t> matches;

            // Extract the target plan's join structure and filters
            std::unordered_map<Subproblem, Subproblem, SubproblemHash> target_structure;
            std::vector<FilterPredicate> target_filters;
            std::unordered_map<Subproblem, std::vector<size_t>, SubproblemHash> target_subproblem_filters;

            // Extract join structure and filters from the plan table
            for (std::size_t i = 1; i < plan_table.size(); ++i)
            {
                Subproblem s(i);
                if (plan_table.has_plan(s))
                {
                    const auto &entry = plan_table[s];

                    // Extract join structure
                    if (s.size() > 1 && entry.left && entry.right)
                    {
                        target_structure[s] = entry.left;
                    }

                    // Extract filter predicates
                    if (entry.data)
                    {
                        if (auto cardinality_data = dynamic_cast<PlanTableEntryCardinalityData *>(entry.data.get()))
                        {
                            if (cardinality_data->data && cardinality_data->data->has_filter)
                            {
                                for (const auto &filter_pred : cardinality_data->data->filter_predicates)
                                {
                                    // Store the filter predicate
                                    size_t pred_idx = target_filters.size();
                                    target_filters.push_back(filter_pred);

                                    // Map subproblem to filter
                                    target_subproblem_filters[s].push_back(pred_idx);
                                }
                            }
                        }
                    }
                }
            }

            // Check each stored plan for a match
            for (size_t i = 0; i < stored_query_plans_.size(); i++)
            {
                const auto &stored_plan = stored_query_plans_[i];

                // First check: Join structure must match
                if (stored_plan.join_structure != target_structure)
                {
                    continue; // Join structure doesn't match
                }

                // Second check: For exact matching, filter predicates must match
                if (exact_match)
                {
                    // Check if filter counts match
                    if (stored_plan.filter_predicates.size() != target_filters.size())
                    {
                        continue;
                    }

                    // Check if filters match exactly (regardless of order)
                    bool filters_match = true;

                    // Check if the same filters apply to the same subproblems
                    for (const auto &[subp, target_indices] : target_subproblem_filters)
                    {
                        auto it = stored_plan.subproblem_filters.find(subp);
                        if (it == stored_plan.subproblem_filters.end() ||
                            it->second.size() != target_indices.size())
                        {
                            filters_match = false;
                            break;
                        }

                        // Check each filter on this subproblem
                        std::vector<FilterPredicate> target_subp_filters;
                        for (size_t idx : target_indices)
                        {
                            target_subp_filters.push_back(target_filters[idx]);
                        }

                        std::vector<FilterPredicate> stored_subp_filters;
                        for (size_t idx : it->second)
                        {
                            stored_subp_filters.push_back(stored_plan.filter_predicates[idx]);
                        }

                        // Compare the filter sets
                        for (const auto &pred : target_subp_filters)
                        {
                            bool found = false;
                            for (const auto &stored_pred : stored_subp_filters)
                            {
                                if (pred == stored_pred)
                                {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found)
                            {
                                filters_match = false;
                                break;
                            }
                        }

                        if (!filters_match)
                            break;
                    }

                    if (!filters_match)
                        continue;
                }

                // If we get here, this plan is a match
                matches.push_back(i);
            }

            if (debug_output_)
            {
                std::cout << "Found " << matches.size() << " matching query plans"
                          << " (exact_match=" << (exact_match ? "true" : "false") << ")" << std::endl;
            }

            return matches;
        }

        /**
         * @brief Get a stored query plan by index
         */
        const StoredQueryPlan *get_stored_query_plan(size_t index) const
        {
            if (index < stored_query_plans_.size())
            {
                return &stored_query_plans_[index];
            }
            return nullptr;
        }

        /**
         * @brief Get the number of stored query plans
         */
        size_t get_stored_query_plan_count() const
        {
            return stored_query_plans_.size();
        }

        /**
         * @brief Look up join cardinality with filter matching
         */
        double lookup_join_cardinality(const Subproblem &left_sp,
                                       const Subproblem &right_sp,
                                       const std::vector<FilterPredicate> *filters,
                                       bool &found) const
        {
            // Simple implementation - look up in our cache
            const Subproblem joined = left_sp | right_sp;

            // Try to find the cardinality for this exact join in our cache
            auto it = subproblem_to_data_.find(joined);
            if (it != subproblem_to_data_.end())
            {
                const auto &data = all_cardinality_data_[it->second];

                // Check if we need to match filters
                if (filters != nullptr)
                {
                    // We need to check if the filters match
                    if (!data->filter_predicates_match(*filters))
                    {
                        found = false;
                        return -1.0;
                    }
                }
                else if (!data->filter_predicates.empty())
                {
                    // There are filters in the stored data but none provided for matching
                    found = false;
                    return -1.0;
                }

                if (data->true_cardinality >= 0)
                {
                    found = true;
                    if (debug_output_)
                    {
                        std::cout << "Found stored cardinality for join: " << data->true_cardinality << std::endl;
                    }
                    return data->true_cardinality;
                }
            }

            found = false;
            return -1.0;
        }

        // Overload that accepts no filters for backward compatibility
        double lookup_join_cardinality(const Subproblem &left_sp,
                                       const Subproblem &right_sp,
                                       bool &found) const
        {
            return lookup_join_cardinality(left_sp, right_sp, nullptr, found);
        }
    };
} // namespace m
