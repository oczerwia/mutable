// TODO
// - Also store the estimated cardinality, so that we can just extract it later
// - Use single method instead of has_stored_cardinality and has_stored_cardinality_range
#pragma once

#include <mutable/IR/PlanTable.hpp>
#include <mutable/IR/Operator.hpp>
#include <mutable/IR/CNF.hpp>
#include <mutable/util/ADT.hpp>
#include <mutable/Options.hpp>


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

    enum class OperatorType
    {
        SCAN,
        FILTER,
        JOIN,
        GROUP_BY,
        PROJECTION,
        LIMIT,
        AGGREGATION,
        OTHER
    };
    inline std::ostream &operator<<(std::ostream &os, const OperatorType &type)
    {
        switch (type)
        {
        case OperatorType::SCAN:
            os << "SCAN";
            break;
        case OperatorType::FILTER:
            os << "FILTER";
            break;
        case OperatorType::JOIN:
            os << "JOIN";
            break;
        case OperatorType::GROUP_BY:
            os << "GROUP_BY";
            break;
        case OperatorType::PROJECTION:
            os << "PROJECTION";
            break;
        case OperatorType::LIMIT:
            os << "LIMIT";
            break;
        case OperatorType::AGGREGATION:
            os << "AGGREGATION";
            break;
        default:
            os << "OTHER";
            break;
        }
        return os;
    }

    /**
     * @brief Complete data structure containing all relevant cardinality information for easy extraction
     */
    struct CardinalityData
    {
        // Cardinality information
        double estimated_cardinality = -1.0;
        double true_cardinality = -1.0;
        std::pair<double, double> estimated_range = {-1.0, -1.0};
        std::size_t operator_order = 0;
        OperatorType operator_type = OperatorType::OTHER;
        std::string operator_name;

        bool has_filter = false;
        bool has_join = false;
        bool has_grouping = false;

        std::set<std::string> filter_strings;
        std::set<std::string> table_names;
        std::set<std::string> group_by_columns;

        // Performance metrics
        double selectivity = -1.0;
        double error_percent = -1.0;

    public:
        double get_cardinality() const { return true_cardinality; }
        bool has_range() const { return estimated_range.first >= 0.0; }
        void set_range(const std::pair<double, double> &range) { estimated_range = range; }
        std::pair<double, double> get_range() const { return estimated_range; }

        bool operator==(const CardinalityData &other) const
        {
            return operator_type == other.operator_type &&
                   table_names == other.table_names &&
                   has_grouping == other.has_grouping &&
                   has_filter == other.has_filter &&
                   filter_strings == other.filter_strings &&
                   group_by_columns == other.group_by_columns;
            // can turn into hash later?
        }

        bool operator!=(const CardinalityData &other) const
        {
            return !(*this == other);
        }
    };


    /**
     * @brief Singleton class that stores cardinality information during query execution
     */
    class CardinalityStorage
    {
    private:

        inline static size_t query_counter_ = 0;

        std::vector<std::shared_ptr<CardinalityData>> stored_cardinalities_;

        bool debug_output_ = true;

        std::vector<std::shared_ptr<CardinalityData>> current_cardinality_data; // This is from the current query, will be reset after each query

        std::vector<std::string> current_table_names;
        std::set<std::string> current_filters_;
        std::set<std::string> current_group_by_columns;

        double cardinality_ = -1.0;
        std::pair<double, double> estimated_range = {-1.0, -1.0};

        bool allow_learning = Options::Get().learn_cardinalities;

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

        /**
         * @brief Get the current table names
         *
         * @return const std::set<std::string>& Reference to current table names
         */
        const std::vector<std::string> &get_current_table_names() const
        {
            return current_table_names;
        }

        /**
         * @brief Set the current GROUP BY columns
         *
         * @param columns The columns used in the GROUP BY clause
         */
        void set_current_group_by_columns(const std::set<std::string> &columns)
        {
            current_group_by_columns = columns;

            if (debug_output_)
            {
                std::cout << "Updated current GROUP BY columns: ";
                for (const auto &col : current_group_by_columns)
                {
                    std::cout << col << " ";
                }
                std::cout << std::endl;
            }
        }

        /**
         * @brief Get the current GROUP BY columns
         *
         * @return const std::set<std::string>& Reference to current GROUP BY columns
         */
        const std::set<std::string> &get_current_group_by_columns() const
        {
            return current_group_by_columns;
        }

        /**
         * @brief Update the current table names set by traversing the plan table
         *
         * @param plan_table The plan table to extract table names from
         */
        void update_current_table_names(const QueryGraph &query_graph)
        {
            std::vector<std::string> table_names;

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
            current_table_names = table_names;

            if (debug_output_)
            {
                std::cout << "Updated current table names from QueryGraph: ";
                for (const auto &name : current_table_names)
                {
                    std::cout << name << " ";
                }
                std::cout << std::endl;
            }
        }

        std::set<std::string> extract_filter_clauses(const cnf::CNF& filter) {
            std::set<std::string> clauses;
            for (const auto& clause : filter) {
                std::ostringstream oss;
                clause.to_sql(oss);
                clauses.insert(oss.str());
            }
            return clauses;
        }

        bool has_stored_cardinality(SmallBitset involved_tables)
        {
            // 1. Iterate over the Bitset
            // 2. get the involved tables from previously determined set of involved tables
            // 3. iterate over the vector of all previously stored queries and search for equal
            //      sets of tables
            if (!allow_learning) { return false; }

            std::set<std::string> table_names;
            for (auto pos : involved_tables)
            {
                if (pos < current_table_names.size())
                {
                    table_names.insert(current_table_names[pos]);
                }
            }
            // could be turned into single if instead of multiple
            for (const auto &stored_cardinality : stored_cardinalities_)
            {
                if (stored_cardinality->table_names == table_names)
                {
                    if (stored_cardinality->filter_strings == current_filters_)
                    {
                        bool group_by_matches =
                            !stored_cardinality->has_grouping ||
                            (stored_cardinality->has_grouping &&
                             stored_cardinality->group_by_columns == current_group_by_columns);

                        if (group_by_matches)
                        {
                            cardinality_ = stored_cardinality->get_cardinality();
                            estimated_range = stored_cardinality->get_range();

                            return true;
                        }
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

        std::pair<double, double> get_stored_cardinality_range() const
        {
            return estimated_range;
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
                    auto clause_strings = extract_filter_clauses(ds.filter());
                    filter_strings.insert(clause_strings.begin(), clause_strings.end());
                }
            }
            current_filters_ = filter_strings;
        }

        struct OperatorCounter {
            static size_t& value() {
                static size_t counter = 0;
                return counter;
            }
        };

        void reset_traverse_counter_() {
            OperatorCounter::value() = 0;
        }


        /** TODO: Check if this is even needed
         * We could use the switch statement to extract the attributes from the operator classes
         * @brief Extract metadata from operators (tables, filters, joins) - SIMPLIFIED
         *
         * @param op The operator to extract metadata from
         * @param data The CardinalityData to update
         */
        void extract_operator_metadata_(const Operator &op, CardinalityData &data)
        {

            if (auto filter_op = dynamic_cast<const FilterOperator *>(&op))
            {
                data.operator_type = OperatorType::FILTER;
                data.has_filter = true;
                auto clause_strings = extract_filter_clauses(filter_op->filter());
                data.filter_strings.insert(clause_strings.begin(), clause_strings.end());
            }
            else if (auto join_op = dynamic_cast<const JoinOperator *>(&op))
            {
                data.operator_type = OperatorType::JOIN;
                data.has_join = true;
            }
            else if (auto group_op = dynamic_cast<const GroupingOperator *>(&op))
            {
                data.operator_type = OperatorType::GROUP_BY;

                std::set<std::string> group_by_cols = extract_group_by_columns(group_op->group_by());
                if (!group_by_cols.empty())
                {
                    data.group_by_columns = group_by_cols;
                    data.has_grouping = true;
                }
            }
            else if (auto scan_op = dynamic_cast<const ScanOperator *>(&op))
            {
                data.operator_type = OperatorType::SCAN;
                // Get table name directly from the scan operator
                const Store &store = scan_op->store();
                const Table &table = store.table();
                data.table_names.insert(std::string(*table.name()));
            }
            else if (dynamic_cast<const ProjectionOperator *>(&op))
            {
                data.operator_type = OperatorType::PROJECTION;
            }
            else if (dynamic_cast<const LimitOperator *>(&op))
            {
                data.operator_type = OperatorType::LIMIT;
            }
            else if (dynamic_cast<const AggregationOperator *>(&op))
            {
                data.operator_type = OperatorType::AGGREGATION;
            }
            else
            {
                data.operator_type = OperatorType::OTHER;
            }

            data.operator_name = typeid(op).name();
        }

        /**
         * @brief Clear current GROUP BY columns (when processing a query without GROUP BY)
         */
        void clear_current_group_by_columns()
        {
            current_group_by_columns.clear();
        }

        /**
         * @brief Extract column names from GROUP BY expressions
         */
        std::set<std::string> extract_group_by_columns(const std::vector<GroupingOperator::group_type> &groups)
        {
            std::set<std::string> group_by_cols;
            for (const auto &[expr, alias] : groups)
            {
                const ast::Expr &expr_ref = expr.get();
                if (auto designator = cast<const ast::Designator>(&expr_ref))
                {
                    std::string col_name;
                    if (designator->has_table_name())
                    {
                        col_name = std::string(*designator->table_name.text) + "." +
                                   std::string(*designator->attr_name.text);
                    }
                    else
                    {
                        col_name = std::string(*designator->attr_name.text);
                    }
                    group_by_cols.insert(col_name);
                }
            }
            return group_by_cols;
        }

        /**
         * @brief Traverses the operator tree DFS and collects cardinality information
         *
         * @param op The current operator being processed
         *
         * DFS is used so that we can propagate filter, group by, ... to parent nodes
         */
        std::shared_ptr<CardinalityData> traverse_operator_tree_(const Operator &op)
        {
            size_t &current_order = OperatorCounter::value();

            // DFS
            std::vector<std::shared_ptr<CardinalityData>> children = {};
            if (auto consumer = dynamic_cast<const Consumer *>(&op))
            {
                for (auto child : consumer->children())
                {
                    auto informed_child = traverse_operator_tree_(*child);
                    children.push_back(informed_child);
                }
            }

            std::shared_ptr<CardinalityData> data = std::make_shared<CardinalityData>();
            data->operator_order = current_order++;

            data->true_cardinality = op.get_emitted_tuples();

            if (op.has_info())
            {
                extract_operator_metadata_(op, *data);
                data->estimated_cardinality = op.info().estimated_cardinality;
                data->estimated_range = op.info().estimated_range;
            }

            if (!children.empty())
            {
                for (const auto child : children)
                {
                    if (child->has_grouping)
                    {
                        data->has_grouping = true;
                        data->group_by_columns.insert(
                            child->group_by_columns.begin(),
                            child->group_by_columns.end());
                    }
                    if (child->has_filter)
                    {
                        data->has_filter = true;
                        data->filter_strings.insert(
                            child->filter_strings.begin(),
                            child->filter_strings.end());
                    }
                    if (child->has_join)
                    {
                        data->has_join = true; // Not really needed
                    }
                    data->table_names.insert(
                        child->table_names.begin(),
                        child->table_names.end());
                }
            }

            current_cardinality_data.push_back(data);
            return data;
        }
        /**
         * @brief Traverse the physical operator tree and collect cardinality information for each operator. 
         * This function is mostly for the debug_print.
         *
         * This function performs a depth-first traversal of the given operator tree (physical plan root),
         * extracting and recording relevant cardinality and metadata for each operator node in the plan.
         * For each operator, it:
         *   - Assigns a unique operator order (DFS pre-order).
         *   - Extracts operator type, estimated and true cardinality, filter and group-by information, and involved tables.
         *   - Propagates filter, group-by, and table information from child operators to parent operators.
         *   - Stores the collected data in the current_cardinality_data vector for later export or learning.
         *
         * This function is typically called once per query, after query execution, and before exporting or learning cardinalities.
         *
         * @param root The root operator of the executed physical plan.
         */
        void map_true_cardinalities_to_logical_plan_(const Operator &root)
        {

            // Clear previous data
            current_cardinality_data.clear();

            reset_traverse_counter_();

            // Traverse and collect cardinalities
            traverse_operator_tree_(root);

            // Debug output
            if (debug_output_)
            {
                std::cout << "\nCollected " << current_cardinality_data.size() << " cardinality entries." << std::endl;
                for (const auto &data : current_cardinality_data)
                {
                    std::cout << "  Subproblem: " << "operator=" << data->operator_type
                              << ", order=" << data->operator_order << ", ";
                    std::cout << "estimated=" << data->estimated_cardinality;
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

                    if (data->has_grouping)
                    {
                        std::cout << "    GROUP BY columns: ";
                        for (const auto &col : data->group_by_columns)
                        {
                            std::cout << col << " ";
                        }
                        std::cout << std::endl;
                    }
                    if (data->has_filter) {
                        std::cout << "    Filter conditions: ";
                        for (const auto &filter : data->filter_strings) {
                            std::cout << filter << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
        }

        void reset_state_for_new_query()
        {
            current_filters_.clear();
            current_group_by_columns.clear();
            current_cardinality_data.clear();
        }

        void clear_stored_operators() {
            current_cardinality_data.clear();
            reset_traverse_counter_();
        }

        /**
         * TODO: Extremely confusing naming conventions, function control flow is hard to grasp
         * @brief Extract true cardinalities and populate reduced query graphs
         *
         * @param root The root operator of the executed physical plan
         */
        void extract_cardinalities_from_execution(const Operator &root)
        {

            query_counter_++; // New select statement

            map_true_cardinalities_to_logical_plan_(root);

            std::vector<CardinalityData> temp_cardinalities;

            for (const auto &cardinality_data : current_cardinality_data)
            {                           
                CardinalityData new_cardinality;


                new_cardinality.table_names.insert(cardinality_data->table_names.begin(),
                                                    cardinality_data->table_names.end());


                for (const auto &filter_str : cardinality_data->filter_strings)
                {
                    new_cardinality.filter_strings.insert(filter_str);
                }
                if (cardinality_data->has_grouping)
                {
                    new_cardinality.group_by_columns = cardinality_data->group_by_columns;
                    new_cardinality.has_grouping = true;
                }

                new_cardinality.operator_type = cardinality_data->operator_type;
                new_cardinality.operator_name = cardinality_data->operator_name;
                new_cardinality.operator_order = cardinality_data->operator_order;

                new_cardinality.true_cardinality = cardinality_data->true_cardinality;
                new_cardinality.set_range(cardinality_data->estimated_range);

                new_cardinality.estimated_cardinality = cardinality_data->estimated_cardinality;

                temp_cardinalities.push_back(new_cardinality);
            }

            for (const CardinalityData &new_cardinality : temp_cardinalities)
            {
                bool found_existing = false;
                for (std::shared_ptr<CardinalityData> &existing_cardinality : stored_cardinalities_)
                {
                    if (existing_cardinality->operator_type == new_cardinality.operator_type &&
                        existing_cardinality->table_names == new_cardinality.table_names &&
                        existing_cardinality->filter_strings == new_cardinality.filter_strings &&
                        existing_cardinality->has_grouping == new_cardinality.has_grouping &&
                        (!existing_cardinality->has_grouping || existing_cardinality->group_by_columns == new_cardinality.group_by_columns))
                    {
                        existing_cardinality->true_cardinality = new_cardinality.true_cardinality;
                        existing_cardinality->estimated_range = new_cardinality.estimated_range;
                        existing_cardinality->estimated_cardinality = new_cardinality.estimated_cardinality;
                        existing_cardinality->group_by_columns = new_cardinality.group_by_columns;
                        found_existing = true;

                        if (debug_output_)
                        {
                            std::cout << "  Updated existing CardinalityData: tables={";
                            for (const auto &table : existing_cardinality->table_names)
                            {
                                std::cout << table << " ";
                            }
                            std::cout << "} Operator type: " << existing_cardinality->operator_type << ", filters={";
                            for (const auto &filter : existing_cardinality->filter_strings)
                            {
                                std::cout << filter << " ";
                            }
                            std::cout << "}";

                            if (existing_cardinality->has_grouping && !existing_cardinality->group_by_columns.empty())
                            {
                                std::cout << ", group_by={";
                                for (const auto &col : existing_cardinality->group_by_columns)
                                {
                                    std::cout << col << " ";
                                }
                                std::cout << "}";
                            }

                            std::cout << ", cardinality=" << existing_cardinality->true_cardinality;
                            if (new_cardinality.has_range())
                            {
                                std::cout << ", range=[" << existing_cardinality->estimated_range.first
                                          << "-" << existing_cardinality->estimated_range.second << "]";
                            }
                            std::cout << std::endl;
                        }
                        break;
                    }
                }

                if (!found_existing && new_cardinality.operator_type != OperatorType::OTHER)
                {
                    stored_cardinalities_.push_back(std::make_shared<CardinalityData>(new_cardinality));

                    if (debug_output_)
                    {
                        std::cout << "  Added new CardinalityData: tables={";
                        for (const auto &table : new_cardinality.table_names)
                        {
                            std::cout << table << " ";
                        }
                        std::cout << "} Operator type: " << new_cardinality.operator_type << ", filters={";
                        for (const auto &filter : new_cardinality.filter_strings)
                        {
                            std::cout << filter << " ";
                        }
                        std::cout << "}";

                        if (new_cardinality.has_grouping && !new_cardinality.group_by_columns.empty())
                        {
                            std::cout << ", group_by={";
                            for (const auto &col : new_cardinality.group_by_columns)
                            {
                                std::cout << col << " ";
                            }
                            std::cout << "}";
                        }

                        std::cout << ", cardinality=" << new_cardinality.true_cardinality;
                        if (new_cardinality.has_range())
                        {
                            std::cout << ", range=[" << new_cardinality.estimated_range.first
                                      << "-" << new_cardinality.estimated_range.second << "]";
                        }
                        std::cout << std::endl;
                    }
                }
            }
        }

        /**
         * TODO: Might be useless, currently used in the group by estimator, better to change later
         * TODO: Can change by storing each operator storage seperatly and each iterator only accesses these instead
         * @brief Check for stored GROUP BY cardinality and update the output model
         *
         * @param G the QueryGraph (needed for compatibility with other methods)
         * @param data_model the input data model containing table information
         * @param groups the GROUP BY expressions
         * @param output_model the output model to update with cardinality
         * @return true if a stored cardinality was found and applied
         */
        bool apply_stored_grouping_cardinality(
            const QueryGraph &G,
            const DataModel &data_model,
            const std::vector<CardinalityEstimator::group_type> &groups,
            DataModel &output_model)
        {
            if (!allow_learning) { return false; }
            // Extract GROUP BY columns using our helper method
            std::set<std::string> group_by_cols = extract_group_by_columns(groups);

            // Set current GROUP BY columns
            set_current_group_by_columns(group_by_cols);

            // Get table names directly from the data model
            const std::set<std::string> &table_names = data_model.original_tables;

            // Check stored cardinality directly using table names
            for (const auto &stored_cardinality : stored_cardinalities_)
            {
                if (stored_cardinality->table_names == table_names &&
                    stored_cardinality->filter_strings == current_filters_)
                {
                    bool group_by_matches =
                        !stored_cardinality->has_grouping ||
                        (stored_cardinality->has_grouping &&
                         stored_cardinality->group_by_columns == current_group_by_columns);

                    if (group_by_matches)
                    {
                        output_model.set_cardinality(stored_cardinality->get_cardinality());
                        if (stored_cardinality->has_range())
                        {
                            output_model.set_range(stored_cardinality->get_range());
                        }
                        return true;
                    }
                }
            }

            return false;
        }

        /**
         * @brief Check for stored AGGREGATION cardinality and update the output model
         *
         * @param G the QueryGraph (needed for compatibility with other methods)
         * @param data_model the input data model containing table information
         * @param output_model the output model to update with cardinality
         * @return true if a stored cardinality was found and applied
         */
        bool apply_stored_aggregation_cardinality(
            const QueryGraph &G,
            const DataModel &data_model,
            DataModel &output_model)
        {
            if (!allow_learning) { return false; }
            // Get table names directly from the data model
            const std::set<std::string> &table_names = data_model.original_tables;
            
            // Check stored cardinality directly using table names
            for (const auto &stored_cardinality : stored_cardinalities_)
            {
                if (stored_cardinality->table_names == table_names &&
                    stored_cardinality->filter_strings == current_filters_ &&
                    stored_cardinality->operator_type == OperatorType::AGGREGATION)
                {
                    // Found a matching entry - update output model
                    output_model.set_cardinality(stored_cardinality->true_cardinality);
                    if (stored_cardinality->has_range()) {
                        output_model.set_range(stored_cardinality->get_range());
                    }
                    
                    if (debug_output_) {
                        std::cout << "  Applied stored aggregation cardinality: " 
                                << stored_cardinality->true_cardinality
                                << " for tables: ";
                        for (const auto &name : table_names) {
                            std::cout << name << " ";
                        }
                        std::cout << std::endl;
                    }
                    
                    return true;
                }
            }
            
            return false;
        }

        /**
         * @brief Check for stored FILTER cardinality and update the output model
         *
         * @param G the QueryGraph (needed for compatibility with other methods)
         * @param data_model the input data model containing table information
         * @param filter the filter condition
         * @param output_model the output model to update with cardinality
         * @return true if a stored cardinality was found and applied
         */
        bool apply_stored_filter_cardinality(
            const QueryGraph &G,
            const DataModel &data_model,
            const cnf::CNF &filter,
            DataModel &output_model)
        {
            if (!allow_learning) { return false; }
            // Early return if filter is empty
            if (filter.empty()) {
                return false;
            }
            
            std::set<std::string> filter_strings = extract_filter_clauses(filter);
            
            // Get table names directly from the data model
            const std::set<std::string> &table_names = data_model.original_tables;
            
            // Check stored cardinality directly using table names and filter
            for (const auto &stored_cardinality : stored_cardinalities_)
            {
                // Match based on table names, operator type, and filter string
                if (stored_cardinality->table_names == table_names &&
                    stored_cardinality->operator_type == OperatorType::FILTER &&
                    stored_cardinality->filter_strings == filter_strings)
                {
                    // Found a matching entry - update output model
                    output_model.set_cardinality(stored_cardinality->true_cardinality);
                    if (stored_cardinality->has_range()) {
                        output_model.set_range(stored_cardinality->get_range());
                    }
                    
                    if (debug_output_) {
                        std::cout << "Stored filters available:" << std::endl;
                        for (const auto& stored : stored_cardinalities_) {
                            if (stored->operator_type == OperatorType::FILTER && 
                                stored->table_names == table_names) {
                                for (const auto& f : stored->filter_strings) {
                                    std::cout << "  '" << f << "'" << std::endl;
                                }
                            }
                        }
                    }
                    
                    return true;
                }
            }
            
            return false;
        }

        /**
         * @brief Export only the current query's cardinality data to CSV
         * 
         * @param filename The CSV file to write to (defaults to "cardinality_data.csv")
         * @return true if export was successful, false otherwise
         */
        inline bool export_to_csv(const std::string& filename = "") 
        {
            if (current_cardinality_data.empty()) {
                return true;
            }

            std::string output_file;
            if (!filename.empty()) {
                output_file = filename;
            } else {
                output_file = std::string("cardinality_data.csv");
            }
            
            std::ofstream csv_file;
            csv_file.open(output_file, std::ios::app);
            
            if (!csv_file.is_open()) {
                std::cerr << "Failed to open CSV file: " << output_file << std::endl;
                return false;
            }
            
            if (csv_file.tellp() == 0) {
                csv_file << "query_id,operator_id,operator_type,tables,est_card,true_card,q_error,filter_conditions,group_by_columns,lower_bound,upper_bound\n";
            }
            
            for (const auto& data : current_cardinality_data) {
                double q_error = -1.0;
                if (data->estimated_cardinality > 0 && data->true_cardinality > 0) {
                    q_error = std::max(data->estimated_cardinality / data->true_cardinality, 
                                    data->true_cardinality / data->estimated_cardinality);
                }
                
                std::string tables = "";
                for (const auto& table : data->table_names) {
                    if (!tables.empty()) tables += "|";
                    tables += table;
                }
                
                std::string filters = "";
                for (const auto& filter : data->filter_strings) {
                    if (!filters.empty()) filters += "|";
                    filters += filter;
                }
                
                std::string group_by = "";
                for (const auto& col : data->group_by_columns) {
                    if (!group_by.empty()) group_by += "|";
                    group_by += col;
                }
                
                csv_file << query_counter_ << ","
                        << data->operator_order << ","
                        << "\"" << data->operator_type << "\","
                        << "\"" << tables << "\","
                        << data->estimated_cardinality << ","
                        << data->true_cardinality << ","
                        << q_error << ","
                        << "\"" << filters << "\","
                        << "\"" << group_by << "\","
                        << data->estimated_range.first << ","
                        << data->estimated_range.second
                        << std::endl;
            }
            
            csv_file.close();
            return true;
        }

    private:
        CardinalityStorage() = default;
    };
} // namespace m
