#include <mutable/catalog/CardinalityEstimator.hpp>

#include "backend/Interpreter.hpp"
#include "catalog/SpnWrapper.hpp"
#include "util/Spn.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <mutable/catalog/Catalog.hpp>
#include <mutable/IR/CNF.hpp>
#include <mutable/IR/Operator.hpp>
#include <mutable/IR/PlanTable.hpp>
#include <mutable/IR/QueryGraph.hpp>
#include <mutable/Options.hpp>
#include <mutable/util/Diagnostic.hpp>
#include <mutable/util/Pool.hpp>
#include <nlohmann/json.hpp>
#include <mutable/catalog/TableStatistics.hpp>
#include <mutable/parse/AST.hpp>

using namespace m;

namespace
{

    namespace options
    {

        std::filesystem::path injected_cardinalities_file;

    }

}

/*======================================================================================================================
 * CardinalityEstimator
 *====================================================================================================================*/

DataModel::~DataModel() {}

std::pair<double, double> DataModel::get_range() const
{
    if (!has_range())
    {
        return {size, size};
    }
    return range;
}

CardinalityEstimator::~CardinalityEstimator() {}

double CardinalityEstimator::predict_number_distinct_values(const DataModel &) const
{
    throw data_model_exception("predicting the number of distinct values is not supported by this data model.");
};

M_LCOV_EXCL_START
void CardinalityEstimator::dump(std::ostream &out) const
{
    print(out);
    out << std::endl;
}

void CardinalityEstimator::dump() const { dump(std::cerr); }
M_LCOV_EXCL_STOP

/*======================================================================================================================
 * CartesianProductEstimator
 *====================================================================================================================*/

std::unique_ptr<DataModel> CartesianProductEstimator::empty_model() const
{
    auto model = std::make_unique<CartesianProductDataModel>();
    model->size = 0;
    return model;
}

std::unique_ptr<DataModel> CartesianProductEstimator::estimate_scan(const QueryGraph &G, Subproblem P) const
{
    M_insist(P.size() == 1, "Subproblem must identify exactly one DataSource");
    auto idx = *P.begin();
    auto &BT = as<const BaseTable>(*G.sources()[idx]);
    auto model = std::make_unique<CartesianProductDataModel>();
    model->size = BT.table().store().num_rows();

    auto stats_ptr = BT.table().statistics();
    model->set_stats(*stats_ptr);
    model->original_tables.insert(stats_ptr->table_name);
    
    return model;
}

std::unique_ptr<DataModel>
CartesianProductEstimator::estimate_filter(const QueryGraph &G, const DataModel &_data, const cnf::CNF &filter) const
{
    /* This model cannot estimate the effects of applying a filter. */
    auto &data = as<const CartesianProductDataModel>(_data);
    auto model = std::make_unique<CartesianProductDataModel>(data);
    
    // Check if we have stored cardinality for this filter+table combination
    if (CardinalityStorage::Get().apply_stored_filter_cardinality(G, data, filter, *model)) {
        return model;
    }

    model->set_stats(data.get_stats());
    model->original_tables = data.original_tables;

    return model;
}

std::unique_ptr<DataModel>
CartesianProductEstimator::estimate_limit(const QueryGraph &, const DataModel &_data, std::size_t limit,
                                          std::size_t offset) const
{
    auto data = as<const CartesianProductDataModel>(_data);
    const std::size_t remaining = offset > data.size ? 0UL : data.size - offset;
    auto model = std::make_unique<CartesianProductDataModel>();
    model->size = std::min(remaining, limit);

    auto stats = model->get_stats();
    stats.row_count = model->size;
    model->set_stats(stats);

    model->original_tables = data.original_tables;

    return model;
}

std::unique_ptr<DataModel>
CartesianProductEstimator::estimate_grouping(const QueryGraph &G, const DataModel &_data,
                                             const std::vector<group_type> &groups) const
{
    auto &data = as<const CartesianProductDataModel>(_data);
    auto model = std::make_unique<CartesianProductDataModel>();

    if (CardinalityStorage::Get().apply_stored_grouping_cardinality(G, data, groups, *model)) {
        return model;
    }
    model->size = data.size;

    model->set_stats(data.get_stats());
    model->original_tables = data.original_tables;
    
    return model;
}

std::unique_ptr<DataModel>
CartesianProductEstimator::estimate_join(const QueryGraph &, const DataModel &_left, const DataModel &_right,
                                         const cnf::CNF &) const
{
    auto left = as<const CartesianProductDataModel>(_left);
    auto right = as<const CartesianProductDataModel>(_right);
    auto model = std::make_unique<CartesianProductDataModel>();

    // Set size as before
    model->size = left.size * right.size;

    // Propagate and merge statistics if available
    auto left_stats = left.get_stats();
    auto right_stats = right.get_stats();
    auto merged_stats = left_stats.merge_for_join(right_stats);
    merged_stats.row_count = model->size;
    model->set_stats(merged_stats);

    // Propagate original tables
    std::set<std::string> all_tables = left.original_tables;
    all_tables.insert(right.original_tables.begin(), right.original_tables.end());
    model->original_tables = all_tables;

    return model;
}

template <typename PlanTable>
std::unique_ptr<DataModel>
CartesianProductEstimator::operator()(estimate_join_all_tag, PlanTable &&PT, const QueryGraph &, Subproblem to_join,
                                      const cnf::CNF &) const
{
    M_insist(not to_join.empty());
    auto model = std::make_unique<CartesianProductDataModel>();
    model->size = 1UL;
    for (auto it = to_join.begin(); it != to_join.end(); ++it)
        model->size *= as<const CartesianProductDataModel>(*PT[it.as_set()].model).size;
    return model;
}

template std::unique_ptr<DataModel>
CartesianProductEstimator::operator()(estimate_join_all_tag, const PlanTableSmallOrDense &, const QueryGraph &,
                                      Subproblem, const cnf::CNF &) const;
template std::unique_ptr<DataModel>
CartesianProductEstimator::operator()(estimate_join_all_tag, const PlanTableLargeAndSparse &, const QueryGraph &,
                                      Subproblem, const cnf::CNF &) const;

std::size_t CartesianProductEstimator::predict_cardinality(const DataModel &data) const
{
    return as<const CartesianProductDataModel>(data).size;
}

M_LCOV_EXCL_START
void CartesianProductEstimator::print(std::ostream &out) const
{
    out << "CartesianProductEstimator - returns size of the Cartesian product of the given subproblems";
}
M_LCOV_EXCL_STOP

/*======================================================================================================================
 * HistogramEstimator
 *====================================================================================================================*/

std::unique_ptr<DataModel> HistogramEstimator::empty_model() const
{
    auto model = std::make_unique<HistogramDataModel>();
    model->size = 0;
    TableStatistics empty_stats;
    empty_stats.row_count = 0;
    model->set_stats(empty_stats);
    return model;
}

std::unique_ptr<DataModel> HistogramEstimator::estimate_scan(const QueryGraph &G, Subproblem P) const
{
    M_insist(P.size() == 1, "Subproblem must identify exactly one DataSource");
    auto idx = *P.begin();
    auto &BT = as<const BaseTable>(*G.sources()[idx]);
    auto model = std::make_unique<HistogramDataModel>();
    model->size = BT.table().store().num_rows();

    // Load table statistics with histograms
    auto stats_ptr = BT.table().statistics();
    model->set_stats(*stats_ptr);
    model->original_tables.insert(stats_ptr->table_name);

    return model;
}

std::unique_ptr<DataModel>
HistogramEstimator::estimate_filter(const QueryGraph &G, const DataModel &_data, const cnf::CNF &filter) const
{
    auto &data = as<const HistogramDataModel>(_data);
    auto result = std::make_unique<HistogramDataModel>(data);

    if (filter.empty())
    {
        return result;
    }

    if (CardinalityStorage::Get().apply_stored_filter_cardinality(G, data, filter, *result)) {
        return result;
    }

    // Apply histogram-based filtering
    auto current_stats = result->get_stats();
    auto filtered_stats = current_stats.filter_by_cnf(filter);
    result->set_stats(filtered_stats);

    result->size = estimate_cardinality_from_histograms(filtered_stats, current_stats.row_count);

    return result;
}

std::unique_ptr<DataModel>
HistogramEstimator::estimate_limit(const QueryGraph &, const DataModel &_data, std::size_t limit,
                                      std::size_t offset) const
{
    auto &data = as<const HistogramDataModel>(_data);
    const std::size_t remaining = offset > data.size ? 0UL : data.size - offset;
    auto result = std::make_unique<HistogramDataModel>(data);
    result->size = std::min(remaining, limit);

    // Update statistics to reflect the limit
    auto stats = result->get_stats();
    stats.row_count = result->size;
    result->set_stats(stats);

    return result;
}

std::unique_ptr<DataModel>
HistogramEstimator::estimate_grouping(const QueryGraph &G, const DataModel &_data,
                                         const std::vector<group_type> &groups) const
{
    auto &data = as<const HistogramDataModel>(_data);
    auto result = std::make_unique<HistogramDataModel>(data);

    if (groups.empty()) {
        result->size = 1;
        return result;
    }

    std::vector<std::string> group_columns;
    for (const auto &[grp, alias] : groups) {
        std::string col_name = extract_column_name_from_expression(grp);
        if (!col_name.empty()) {
            group_columns.push_back(col_name);
        }
    }

    auto current_stats = result->get_stats();
    result->size = current_stats.estimate_group_by_cardinality(group_columns);
    
    auto grouped_stats = current_stats.apply_group_by(group_columns);
    
    grouped_stats = grouped_stats.rescale_histograms_to_cardinality(result->size);
    
    result->set_stats(grouped_stats);
    
    return result;
}


std::unique_ptr<DataModel>
HistogramEstimator::estimate_join(const QueryGraph &G, const DataModel &_left, const DataModel &_right,
                                     const cnf::CNF &condition) const
{
    auto &left = as<const HistogramDataModel>(_left);
    auto &right = as<const HistogramDataModel>(_right);
    auto result = std::make_unique<HistogramDataModel>();

    auto left_stats = left.get_stats();
    auto right_stats = right.get_stats();

    std::vector<std::pair<std::string, std::string>> join_pairs;
    for (const auto &join : G.joins()) {
        auto join_condition = join->condition();
        if (!join_condition.empty()) {
            auto join_columns = join_condition.get_join_columns();
            for (const auto &[table, columns] : join_columns) {
                if (!columns.empty()) {
                    for (const auto &[other_table, other_columns] : join_columns) {
                        if (table != other_table && !other_columns.empty()) {
                            std::string left_col = table + "." + *columns.begin();
                            std::string right_col = other_table + "." + *other_columns.begin();
                            join_pairs.emplace_back(left_col, right_col);
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    
    if (join_pairs.empty()) {
        result->size = left.size * right.size;
    } else {
        bool found_valid_histogram = false;
        
        for (const auto &[left_col, right_col] : join_pairs) {
            auto join_hist = left_stats.multiply_histograms(left_col, right_col);
            if (join_hist.is_valid()) {
                result->size = std::accumulate(join_hist.bins.begin(), join_hist.bins.end(), 0UL);
                found_valid_histogram = true;
                break;
            }
        }
        
        if (!found_valid_histogram) {
            result->size = left.size * right.size;
        }
    }

    auto merged = left_stats.merge_for_join(right_stats);
    merged.row_count = result->size;
    
    merged = merged.rescale_histograms_to_cardinality(result->size);
    
    std::set<std::string> all_tables = left.original_tables;
    all_tables.insert(right.original_tables.begin(), right.original_tables.end());
    result->original_tables = all_tables;
    
    result->set_stats(merged);
    return result; // there is no joined column on what we can run the next join in the result section
    // NOTE: We will keep it like that and write it down as an assumption
}

template <typename PlanTable>
std::unique_ptr<DataModel>
HistogramEstimator::operator()(estimate_join_all_tag, PlanTable &&PT, const QueryGraph &G, Subproblem to_join,
                                  const cnf::CNF &condition) const
{
    M_insist(not to_join.empty());

    if (to_join.size() == 1)
    {
        auto &single_model = as<const HistogramDataModel>(*PT[to_join.begin().as_set()].model);
        return std::make_unique<HistogramDataModel>(single_model);
    }

    // Join multiple tables iteratively
    auto it = to_join.begin();
    auto &first_model = as<const HistogramDataModel>(*PT[it.as_set()].model);
    auto result_model = std::make_unique<HistogramDataModel>(first_model);
    ++it;

    for (; it != to_join.end(); ++it)
    {
        // Fix: Use reference correctly
        auto &right_model = as<const HistogramDataModel>(*PT[it.as_set()].model);
        auto joined = estimate_join(G, *result_model, right_model, condition);
        result_model = std::unique_ptr<HistogramDataModel>(static_cast<HistogramDataModel *>(joined.release()));
    }

    return result_model;
}

template std::unique_ptr<DataModel>
HistogramEstimator::operator()(estimate_join_all_tag, const PlanTableSmallOrDense &, const QueryGraph &,
                                  Subproblem, const cnf::CNF &) const;
template std::unique_ptr<DataModel>
HistogramEstimator::operator()(estimate_join_all_tag, const PlanTableLargeAndSparse &, const QueryGraph &,
                                  Subproblem, const cnf::CNF &) const;

std::size_t HistogramEstimator::predict_cardinality(const DataModel &data) const
{
    return as<const HistogramDataModel>(data).size;
}

// Helper method implementations for HistogramEstimator
std::size_t HistogramEstimator::estimate_cardinality_from_histograms(const TableStatistics &filtered_stats, std::size_t original_size) const
{
    if (filtered_stats.histograms.empty())
    {
        return 1;
    }

    // TODO: 
    // Via independence assumption, we should adjust every histogram after a filter
    // by rescaling the histogram to the new number of rows and cutting the NDV appropriatly 
    // (assume equal distribution of NDV per bucket)
    // Then we should just return the size after all filters (which should be equal on all histograms)
    // This means that this method should just return the size of any histogram within the table statistics

    std::size_t min_cardinality = std::numeric_limits<std::size_t>::max();
    bool found_valid = false;

    for (const auto &[col_name, hist] : filtered_stats.histograms)
    {
        if (hist.is_valid() && hist.total_count > 0)
        {
            min_cardinality = std::min(min_cardinality, hist.total_count);
            found_valid = true;
        }
    }

    return found_valid ? min_cardinality : original_size;
}

std::vector<std::pair<std::string, std::string>> HistogramEstimator::extract_equi_join_columns(const cnf::CNF &condition) const
{
    std::vector<std::pair<std::string, std::string>> join_pairs;

    for (const auto &clause : condition)
    {
        if (clause.size() != 1)
            continue; // Skip complex clauses

        const auto &predicate = clause[0];
        if (predicate.negative())
            continue; // Skip negated predicates

        if (auto binary_expr = cast<const ast::BinaryExpr>(&predicate.expr()))
        {
            if (binary_expr->op().type == TK_EQUAL)
            {
                // Check if both sides are designators (table.column)
                auto left_des = cast<const ast::Designator>(binary_expr->lhs.get());
                auto right_des = cast<const ast::Designator>(binary_expr->rhs.get());

                if (left_des && right_des &&
                    left_des->has_table_name() && right_des->has_table_name())
                {

                    std::string left_col = std::string(*left_des->table_name.text) + "." +
                                           std::string(*left_des->attr_name.text);
                    std::string right_col = std::string(*right_des->table_name.text) + "." +
                                            std::string(*right_des->attr_name.text);

                    join_pairs.emplace_back(left_col, right_col);
                }
            }
        }
    }

    return join_pairs;
}

std::string HistogramEstimator::extract_column_name_from_expression(const ast::Expr &expr) const
{
    if (auto designator = cast<const ast::Designator>(&expr)) {
        if (designator->has_table_name()) {
            return std::string(*designator->table_name.text) + "." + *designator->attr_name.text;
        } else {
            return std::string(*designator->attr_name.text);
        }
    }
    return ""; // Add fallback for other expression types
}


M_LCOV_EXCL_START
void HistogramEstimator::print(std::ostream &out) const
{
    out << "HistogramEstimator - returns size based on the histogram of the given subproblems";
}
M_LCOV_EXCL_STOP

/*======================================================================================================================
 * RangeCartesianProductEstimator
 *====================================================================================================================*/

std::unique_ptr<DataModel> RangeCartesianProductEstimator::empty_model() const
{
    auto model = std::make_unique<RangeCartesianProductDataModel>();
    model->set_cardinality(0);
    model->set_range({0.0, 0.0});
    return model;
}

std::unique_ptr<DataModel> RangeCartesianProductEstimator::estimate_scan(const QueryGraph &G, Subproblem P) const
{
    M_insist(P.size() == 1, "Subproblem must identify exactly one DataSource");
    auto idx = *P.begin();
    auto &BT = as<const BaseTable>(*G.sources()[idx]);
    double num_rows = static_cast<double>(BT.table().store().num_rows());

    auto model = std::make_unique<RangeCartesianProductDataModel>(num_rows);
    model->set_range({num_rows, num_rows});
    return model;
}

std::unique_ptr<DataModel>
RangeCartesianProductEstimator::estimate_filter(const QueryGraph &G, const DataModel &data, const cnf::CNF &filter) const
{
    auto &model = as<const RangeCartesianProductDataModel>(data);
    auto result = std::make_unique<RangeCartesianProductDataModel>();

    if (filter.empty())
    {
        result->set_cardinality(model.size);
        if (model.has_range())
            result->set_range(model.get_range());
        return result;
    }

    // TODO: Change this to some naive min max range estimation or just return the full result
    auto range = model.has_range() ? model.get_range() : std::make_pair(double(model.size), double(model.size));

    // Simple filter model: reduce by 30-70% (add range uncertainty)
    double low_selectivity = 0.3;  // reduce by at most 70%
    double high_selectivity = 0.7; // reduce by at least 30%

    result->set_cardinality(model.size * high_selectivity);
    result->set_range({range.first * low_selectivity, range.second * high_selectivity});
    return result;
}

std::unique_ptr<DataModel>
RangeCartesianProductEstimator::estimate_limit(const QueryGraph &G, const DataModel &data, std::size_t limit,
                                               std::size_t offset) const
{
    auto &model = as<const RangeCartesianProductDataModel>(data);
    auto result = std::make_unique<RangeCartesianProductDataModel>();

    std::size_t size = std::min(limit, model.size);
    if (offset < model.size)
        size = std::min(size, model.size - offset);
    else
        size = 0;

    result->set_cardinality(size);

    if (model.has_range())
    {
        auto range = model.get_range();
        double min_size = std::min(double(limit), std::max(0.0, range.first - offset));
        double max_size = std::min(double(limit), std::max(0.0, range.second - offset));
        result->set_range({min_size, max_size});
    }
    else
    {
        result->set_range({double(size), double(size)});
    }
    return result;
}

std::unique_ptr<DataModel>
RangeCartesianProductEstimator::estimate_grouping(const QueryGraph &G, const DataModel &data,
                                                  const std::vector<group_type> &groups) const
{
    auto &model = as<const RangeCartesianProductDataModel>(data);
    auto result = std::make_unique<RangeCartesianProductDataModel>();


    std::size_t size = std::min(model.size, std::size_t(groups.empty() ? 1 : model.size));
    result->set_cardinality(size);

    if (model.has_range())
    {
        auto range = model.get_range();
        // For grouping, range depends on group key selectivity
        // Worst case: all values are unique (upper bound = input size)
        // Best case: all values are the same (lower bound = 1)
        double min_size = groups.empty() ? 1.0 : 1.0;
        double max_size = groups.empty() ? 1.0 : range.second;
        result->set_range({min_size, max_size});
    }
    else
    {
        result->set_range({double(size), double(size)});
    }

    if (groups.empty())
    {
        result->set_range({double(1.0), double(1.0)});
    }
    return result;
}

std::unique_ptr<DataModel>
RangeCartesianProductEstimator::estimate_join(const QueryGraph &G, const DataModel &left, const DataModel &right,
                                              const cnf::CNF &condition) const
{
    auto &left_model = as<const RangeCartesianProductDataModel>(left);
    auto &right_model = as<const RangeCartesianProductDataModel>(right);
    auto result = std::make_unique<RangeCartesianProductDataModel>();

    // Get point estimates
    result->set_cardinality(left_model.size * right_model.size);

    // Get ranges
    std::pair<double, double> left_range = left_model.has_range() ? left_model.get_range() : std::make_pair(double(left_model.size), double(left_model.size));
    std::pair<double, double> right_range = right_model.has_range() ? right_model.get_range() : std::make_pair(double(right_model.size), double(right_model.size));

    if (condition.empty())
    {
        // Cartesian product with ranges
        double min_result = left_range.first * right_range.first;
        double max_result = left_range.second * right_range.second;
        result->set_range({min_result, max_result});
    }
    else
    {
        // TODO: Need to remove this first
        double low_selectivity = 0.001; // Highly selective join possible (0.1%)
        double high_selectivity = 0.1;  // Less selective join possible (10%)

        double min_result = std::min(left_range.first, right_range.first) * low_selectivity;
        double max_result = left_range.second * right_range.second * high_selectivity;
        result->set_range({min_result, max_result});
    }

    return result;
}

template <typename PlanTable>
std::unique_ptr<DataModel>
RangeCartesianProductEstimator::operator()(estimate_join_all_tag, PlanTable &&PT, const QueryGraph &G, Subproblem to_join,
                                           const cnf::CNF &condition) const
{
    M_insist(not to_join.empty());
    auto model = std::make_unique<RangeCartesianProductDataModel>();
    model->set_cardinality(1UL);
    model->set_range({1.0, 1.0});

    std::pair<double, double> total_range = {1.0, 1.0};

    for (auto it = to_join.begin(); it != to_join.end(); ++it)
    {
        auto submodel = as<const RangeCartesianProductDataModel>(*PT[it.as_set()].model);
        model->set_cardinality(model->size * submodel.size);

        // Multiply ranges
        if (submodel.has_range())
        {
            auto range = submodel.get_range();
            total_range.first *= range.first;
            total_range.second *= range.second;
        }
        else
        {
            total_range.first *= submodel.size;
            total_range.second *= submodel.size;
        }
    }

    model->set_range(total_range);
    return model;
}

std::size_t RangeCartesianProductEstimator::predict_cardinality(const DataModel &data) const
{
    auto &model = as<const RangeCartesianProductDataModel>(data);
    // Use upper bound for conservative estimate
    if (model.has_range())
        return model.get_range().second;
    return model.size;
}

void RangeCartesianProductEstimator::print(std::ostream &out) const
{
    out << "RangeCartesianProductEstimator";
}

/*======================================================================================================================
 * InjectionCardinalityEstimator
 *====================================================================================================================*/

/*----- Constructors -------------------------------------------------------------------------------------------------*/

InjectionCardinalityEstimator::InjectionCardinalityEstimator(ThreadSafePooledString name_of_database)
    : fallback_(name_of_database)
{
    Diagnostic diag(Options::Get().has_color, std::cout, std::cerr);
    Position pos("InjectionCardinalityEstimator");

    if (options::injected_cardinalities_file.empty())
    {
        std::cout << "No injection file was passed.\n";
    }
    else
    {
        std::ifstream in(options::injected_cardinalities_file);
        if (in)
        {
            read_json(diag, in, name_of_database);
        }
        else
        {
            diag.w(pos) << "Could not open file " << options::injected_cardinalities_file << ".\n"
                        << "A dummy estimator will be used to do estimations.\n";
        }
    }
}

InjectionCardinalityEstimator::InjectionCardinalityEstimator(Diagnostic &diag, ThreadSafePooledString name_of_database,
                                                             std::istream &in)
    : fallback_(name_of_database)
{
    read_json(diag, in, name_of_database);
}

void InjectionCardinalityEstimator::read_json(Diagnostic &diag, std::istream &in,
                                              const ThreadSafePooledString &name_of_database)
{
    Catalog &C = Catalog::Get();
    Position pos("InjectionCardinalityEstimator");
    std::string prev_relation;

    using json = nlohmann::json;
    json cardinalities;
    try
    {
        in >> cardinalities;
    }
    catch (json::parse_error parse_error)
    {
        diag.w(pos) << "The file could not be parsed as json. Parser error output:\n"
                    << parse_error.what() << "\n"
                    << "A dummy estimator will be used to do estimations.\n";
        return;
    }
    json *database_entry;
    try
    {
        database_entry = &cardinalities.at(*name_of_database); // throws if key does not exist
    }
    catch (json::out_of_range &out_of_range)
    {
        diag.w(pos) << "No entry for the db " << name_of_database << " in the file.\n"
                    << "A dummy estimator will be used to do estimations.\n";
        return;
    }
    cardinality_table_.reserve(database_entry->size());
    std::vector<std::string> names;
    for (auto &subproblem_entry : *database_entry)
    {
        json *relations_array;
        json *size;
        try
        {
            relations_array = &subproblem_entry.at("relations");
            size = &subproblem_entry.at("size");
        }
        catch (json::exception &exception)
        {
            diag.w(pos) << "The entry " << subproblem_entry << " for the db \"" << name_of_database << "\""
                        << " does not have the required form of {\"relations\": ..., \"size\": ... } "
                        << "and will thus be ignored.\n";
            continue;
        }

        names.clear();
        for (auto it = relations_array->begin(); it != relations_array->end(); ++it)
            names.emplace_back(it->get<std::string>());
        std::sort(names.begin(), names.end());

        buf_.clear();
        for (auto it = names.begin(); it != names.end(); ++it)
        {
            if (it != names.begin())
                buf_.emplace_back('$');
            buf_append(*it);
        }
        buf_.emplace_back(0);
        ThreadSafePooledString str = C.pool(buf_view());
        auto res = cardinality_table_.emplace(std::move(str), *size);
        M_insist(res.second, "insertion must not fail as we do not allow for duplicates in the input file");
    }
}

/*----- Model calculation --------------------------------------------------------------------------------------------*/

std::unique_ptr<DataModel> InjectionCardinalityEstimator::empty_model() const
{
    return std::make_unique<InjectionCardinalityDataModel>(Subproblem(), 0);
}

std::unique_ptr<DataModel> InjectionCardinalityEstimator::estimate_scan(const QueryGraph &G, Subproblem P) const
{
    M_insist(P.size() == 1);
    const auto idx = *P.begin();
    auto &DS = *G.sources()[idx];

    if (auto it = cardinality_table_.find(DS.name().assert_not_none()); it != cardinality_table_.end())
    {
        return std::make_unique<InjectionCardinalityDataModel>(P, it->second);
    }
    else
    {
        /* no match, fall back */
        auto fallback_model = fallback_.estimate_scan(G, P);
        return std::make_unique<InjectionCardinalityDataModel>(P, fallback_.predict_cardinality(*fallback_model));
    }
}

std::unique_ptr<DataModel>
InjectionCardinalityEstimator::estimate_filter(const QueryGraph &, const DataModel &_data, const cnf::CNF &) const
{
    /* This model cannot estimate the effects of applying a filter. */
    auto &data = as<const InjectionCardinalityDataModel>(_data);
    return std::make_unique<InjectionCardinalityDataModel>(data); // copy
}

std::unique_ptr<DataModel>
InjectionCardinalityEstimator::estimate_limit(const QueryGraph &, const DataModel &_data, std::size_t limit,
                                              std::size_t offset) const
{
    auto &data = as<const InjectionCardinalityDataModel>(_data);
    const std::size_t remaining = offset > data.size_ ? 0UL : data.size_ - offset;
    return std::make_unique<InjectionCardinalityDataModel>(data.subproblem_, std::min(remaining, limit));
}

std::unique_ptr<DataModel>
InjectionCardinalityEstimator::estimate_grouping(const QueryGraph &, const DataModel &_data,
                                                 const std::vector<group_type> &exprs) const
{
    auto &data = as<const InjectionCardinalityDataModel>(_data);

    if (exprs.empty())
        return std::make_unique<InjectionCardinalityDataModel>(data.subproblem_, 1); // single group

    /* Combine grouping keys into an identifier. */
    oss_.str("");
    oss_ << "g";
    for (auto [grp, alias] : exprs)
    {
        oss_ << '#';
        if (alias.has_value())
            oss_ << alias;
        else
            oss_ << grp.get();
    }
    ThreadSafePooledString id = Catalog::Get().pool(oss_.str().c_str());

    if (auto it = cardinality_table_.find(id); it != cardinality_table_.end())
    {
        /* Clamp injected cardinality to at most the cardinality of the grouping's child since it cannot produce more
         * tuples than it receives. */
        return std::make_unique<InjectionCardinalityDataModel>(data.subproblem_, std::min(it->second, data.size_));
    }
    else
    {
        /* This model cannot estimate the effects of grouping. */
        return std::make_unique<InjectionCardinalityDataModel>(data); // copy
    }
}

std::unique_ptr<DataModel>
InjectionCardinalityEstimator::estimate_join(const QueryGraph &G, const DataModel &_left, const DataModel &_right,
                                             const cnf::CNF &condition) const
{
    auto &left = as<const InjectionCardinalityDataModel>(_left);
    auto &right = as<const InjectionCardinalityDataModel>(_right);

    const Subproblem subproblem = left.subproblem_ | right.subproblem_;
    ThreadSafePooledString id = make_identifier(G, subproblem);

    /* Lookup cardinality in table. */
    if (auto it = cardinality_table_.find(id); it != cardinality_table_.end())
    {
        /* Clamp injected cardinality to at most the cardinality of the cartesian product of the join's children
         * since it cannot produce more tuples than that. */
        const std::size_t max_cardinality = left.size_ * right.size_;
        return std::make_unique<InjectionCardinalityDataModel>(subproblem, std::min(it->second, max_cardinality));
    }
    else
    {
        /* Fallback to CartesianProductEstimator. */
        if (not Options::Get().quiet)
            std::cerr << "warning: failed to estimate the join of " << left.subproblem_ << " and " << right.subproblem_
                      << '\n';
        auto left_fallback = std::make_unique<CartesianProductEstimator::CartesianProductDataModel>();
        left_fallback->size = left.size_;
        auto right_fallback = std::make_unique<CartesianProductEstimator::CartesianProductDataModel>();
        right_fallback->size = right.size_;
        auto fallback_model = fallback_.estimate_join(G, *left_fallback, *right_fallback, condition);
        return std::make_unique<InjectionCardinalityDataModel>(subproblem,
                                                               fallback_.predict_cardinality(*fallback_model));
    }
}

template <typename PlanTable>
std::unique_ptr<DataModel>
InjectionCardinalityEstimator::operator()(estimate_join_all_tag, PlanTable &&PT, const QueryGraph &G,
                                          Subproblem to_join, const cnf::CNF &) const
{
    ThreadSafePooledString id = make_identifier(G, to_join);
    if (auto it = cardinality_table_.find(id); it != cardinality_table_.end())
    {
        /* Clamp injected cardinality to at most the cardinality of the cartesian product of the join's children
         * since it cannot produce more tuples than that. */
        std::size_t max_cardinality = 1;
        for (auto it = to_join.begin(); it != to_join.end(); ++it)
            max_cardinality *= as<const InjectionCardinalityDataModel>(*PT[it.as_set()].model).size_;
        return std::make_unique<InjectionCardinalityDataModel>(to_join, std::min(it->second, max_cardinality));
    }
    else
    {
        /* Fallback to cartesian product. */
        if (not Options::Get().quiet)
            std::cerr << "warning: failed to estimate the join of all data sources in " << to_join << '\n';
        auto ds_it = to_join.begin();
        std::size_t size = as<const InjectionCardinalityDataModel>(*PT[ds_it.as_set()].model).size_;
        for (; ds_it != to_join.end(); ++ds_it)
            size *= as<const InjectionCardinalityDataModel>(*PT[ds_it.as_set()].model).size_;
        return std::make_unique<InjectionCardinalityDataModel>(to_join, size);
    }
}

template std::unique_ptr<DataModel>
InjectionCardinalityEstimator::operator()(estimate_join_all_tag, const PlanTableSmallOrDense &, const QueryGraph &,
                                          Subproblem, const cnf::CNF &) const;

template std::unique_ptr<DataModel>
InjectionCardinalityEstimator::operator()(estimate_join_all_tag, const PlanTableLargeAndSparse &, const QueryGraph &,
                                          Subproblem, const cnf::CNF &) const;

std::size_t InjectionCardinalityEstimator::predict_cardinality(const DataModel &data) const
{
    return as<const InjectionCardinalityDataModel>(data).size_;
}

M_LCOV_EXCL_START
void InjectionCardinalityEstimator::print(std::ostream &out) const
{
    constexpr uint32_t max_rows_printed = 100; /// Number of rows of the cardinality_table printed
    std::size_t sub_len = 13;                  /// Length of Subproblem column
    for (auto &entry : cardinality_table_)
        sub_len = std::max(sub_len, strlen(*entry.first));

    out << std::left << "InjectionCardinalityEstimator\n"
        << std::setw(sub_len) << "Subproblem" << "Size" << "\n"
        << std::right;

    /* ------- Print maximum max_rows_printed rows of the cardinality_table_ */
    uint32_t counter = 0;
    for (auto &entry : cardinality_table_)
    {
        if (counter >= max_rows_printed)
            break;
        out << std::left << std::setw(sub_len) << entry.first << entry.second << "\n";
        counter++;
    }
}
M_LCOV_EXCL_STOP

ThreadSafePooledString InjectionCardinalityEstimator::make_identifier(const QueryGraph &G, const Subproblem S) const
{
    auto &C = Catalog::Get();
    static thread_local std::vector<ThreadSafePooledString> names;
    names.clear();
    for (auto id : S)
        names.emplace_back(G.sources()[id]->name());
    std::sort(names.begin(), names.end(), [](auto lhs, auto rhs)
              { return strcmp(*lhs, *rhs) < 0; });

    buf_.clear();
    for (auto it = names.begin(); it != names.end(); ++it)
    {
        if (it != names.begin())
            buf_.emplace_back('$');
        buf_append(**it);
    }

    buf_.emplace_back(0);
    return C.pool(buf_view());
}

/*======================================================================================================================
 * SpnEstimator
 *====================================================================================================================*/

namespace
{

    /** Visitor to translate a CNF to an Spn filter. Only consider sargable expressions. */
    struct FilterTranslator : ast::ConstASTExprVisitor
    {
        Catalog &C = Catalog::Get();
        ThreadSafePooledString attribute;
        float value;
        Spn::SpnOperator op;

        FilterTranslator() : attribute(C.pool("")), value(0), op(Spn::EQUAL) {}

        using ConstASTExprVisitor::operator();

        void operator()(const ast::Designator &designator) { attribute = designator.attr_name.text.assert_not_none(); }

        void operator()(const ast::Constant &constant)
        {
            auto val = Interpreter::eval(constant);

            visit(overloaded{
                      [&val, this](const Numeric &numeric)
                      {
                          switch (numeric.kind)
                          {
                          case Numeric::N_Int:
                              value = float(val.as_i());
                              break;
                          case Numeric::N_Float:
                              value = val.as_f();
                              break;
                          case Numeric::N_Decimal:
                              value = float(val.as_d());
                              break;
                          }
                      },
                      [this](const NoneType &)
                      { op = Spn::IS_NULL; },
                      [](auto &&)
                      { M_unreachable("Unsupported type."); },
                  },
                  *constant.type());
        }

        void operator()(const ast::BinaryExpr &binary_expr)
        {
            switch (binary_expr.op().type)
            {
            case TK_EQUAL:
                op = Spn::EQUAL;
                break;
            case TK_LESS:
                op = Spn::LESS;
                break;
            case TK_LESS_EQUAL:
                op = Spn::LESS_EQUAL;
                break;
            case TK_GREATER:
                op = Spn::GREATER;
                break;
            case TK_GREATER_EQUAL:
                op = Spn::GREATER_EQUAL;
                break;
            default:
                M_unreachable("Operator can't be handled");
            }

            (*this)(*binary_expr.lhs);
            (*this)(*binary_expr.rhs);
        }

        void operator()(const ast::ErrorExpr &) { /* nothing to be done */ }
        void operator()(const ast::FnApplicationExpr &) { /* nothing to be done */ }
        void operator()(const ast::UnaryExpr &) { /* nothing to be done */ }
        void operator()(const ast::QueryExpr &) { /* nothing to be done */ }
    };

    /** Visitor to translate a CNF join condition (get the two identifiers of an equi Join). */
    struct JoinTranslator : ast::ConstASTExprVisitor
    {

        std::vector<std::pair<ThreadSafePooledString, ThreadSafePooledString>> join_designator;

        using ConstASTExprVisitor::operator();

        void operator()(const ast::Designator &designator)
        {
            join_designator.emplace_back(designator.table_name.text, designator.attr_name.text);
        }

        void operator()(const ast::BinaryExpr &binary_expr)
        {

            if (binary_expr.op().type != m::TK_EQUAL)
            {
                M_unreachable("Operator can't be handled");
            }

            (*this)(*binary_expr.lhs);
            (*this)(*binary_expr.rhs);
        }

        void operator()(const ast::Constant &) { /* nothing to be done */ }
        void operator()(const ast::ErrorExpr &) { /* nothing to be done */ }
        void operator()(const ast::FnApplicationExpr &) { /* nothing to be done */ }
        void operator()(const ast::UnaryExpr &) { /* nothing to be done */ }
        void operator()(const ast::QueryExpr &) { /* nothing to be done */ }
    };

}

SpnEstimator::~SpnEstimator()
{
    for (auto &e : table_to_spn_)
        delete e.second;
}

void SpnEstimator::learn_spns() { table_to_spn_ = SpnWrapper::learn_spn_database(name_of_database_); }

void SpnEstimator::learn_new_spn(const ThreadSafePooledString &name_of_table)
{
    table_to_spn_.emplace(
        name_of_table,
        new SpnWrapper(SpnWrapper::learn_spn_table(name_of_database_, name_of_table)));
}

std::pair<unsigned, bool> SpnEstimator::find_spn_id(const SpnDataModel &data, SpnJoin &join)
{
    /* we only have a single spn */
    ThreadSafePooledString table_name = data.spns_.begin()->first;
    auto &attr_to_id = data.spns_.begin()->second.get().get_attribute_to_id();

    unsigned spn_id = 0;
    bool is_primary_key = false;

    /* check which identifier of the join corresponds to the table of the spn, get the corresponding attribute id */
    auto find_iter = table_name == join.first.first ? attr_to_id.find(join.first.second) : attr_to_id.find(join.second.second);

    if (find_iter != attr_to_id.end())
    {
        spn_id = find_iter->second;
    }
    else
    {
        is_primary_key = true;
    }

    return {spn_id, is_primary_key};
}

std::size_t SpnEstimator::max_frequency(const SpnDataModel &data, SpnJoin &join)
{
    auto [spn_id, is_primary_key] = find_spn_id(data, join);

    /* maximum frequency is only computed on data models which only have one Spn */
    const SpnWrapper &spn = data.spns_.begin()->second.get();

    return is_primary_key ? 1 : data.num_rows_ / spn.estimate_number_distinct_values(spn_id);
}

std::size_t SpnEstimator::max_frequency(const SpnDataModel &data, const ThreadSafePooledString &attribute)
{
    /* maximum frequency is only computed on data models which only have one Spn */
    const SpnWrapper &spn = data.spns_.begin()->second.get();
    auto &attr_to_id = spn.get_attribute_to_id();
    auto find_iter = attr_to_id.find(attribute);

    return find_iter == attr_to_id.end() ? 1 : data.num_rows_ / spn.estimate_number_distinct_values(find_iter->second);
}

/*----- Model calculation --------------------------------------------------------------------------------------------*/

std::unique_ptr<DataModel> SpnEstimator::empty_model() const
{
    return std::make_unique<SpnDataModel>(table_spn_map(), 0);
}

std::unique_ptr<DataModel> SpnEstimator::estimate_scan(const QueryGraph &G, Subproblem P) const
{
    M_insist(P.size() == 1);
    const auto idx = *P.begin();
    auto &BT = as<const BaseTable>(*G.sources()[idx]);
    /* get the Spn corresponding for the table to scan */
    if (auto it = table_to_spn_.find(BT.name().assert_not_none()); it != table_to_spn_.end())
    {
        table_spn_map spns;
        const SpnWrapper &spn = *it->second;
        spns.emplace(BT.name(), spn);
        return std::make_unique<SpnDataModel>(std::move(spns), spn.num_rows());
    }
    else
    {
        throw data_model_exception("Table does not exist.");
    }
}

std::unique_ptr<DataModel>
SpnEstimator::estimate_filter(const QueryGraph &, const DataModel &_data, const cnf::CNF &filter) const
{
    auto &data = as<const SpnDataModel>(_data);
    M_insist(data.spns_.size() == 1);
    auto new_data = std::make_unique<SpnDataModel>(data);
    auto &spn = new_data->spns_.begin()->second.get();
    auto &attribute_to_id = spn.get_attribute_to_id();

    Spn::Filter translated_filter;
    FilterTranslator ft;

    /* only consider clauses with one element, since Spns cannot estimate disjunctions */
    for (auto &clause : filter)
    {
        M_insist(clause.size() == 1);
        ft(*clause[0]);
        unsigned spn_id;

        if (auto it = attribute_to_id.find(ft.attribute); it != attribute_to_id.end())
        {
            spn_id = it->second;
        }
        else
        {
            throw data_model_exception("Attribute does not exist.");
        }

        translated_filter.emplace(spn_id, std::make_pair(ft.op, ft.value));
    }

    /* Save number of rows in the newly constructed data model with the filter applied */
    new_data->num_rows_ = float(new_data->num_rows_) * spn.likelihood(translated_filter);
    return new_data;
}

std::unique_ptr<DataModel>
SpnEstimator::estimate_limit(const QueryGraph &, const DataModel &data, std::size_t limit, std::size_t) const
{
    auto model = std::make_unique<SpnDataModel>(as<const SpnDataModel>(data));
    model->num_rows_ = std::min(model->num_rows_, limit);
    return model;
}

std::unique_ptr<DataModel> SpnEstimator::estimate_grouping(const QueryGraph &, const DataModel &data,
                                                           const std::vector<group_type> &groups) const
{
    auto model = std::make_unique<SpnDataModel>(as<const SpnDataModel>(data));
    std::size_t num_rows = 1;
    for (auto [grp, alias] : groups)
    {
        auto designator = cast<const ast::Designator>(&grp.get());
        if (not designator)
            throw data_model_exception("SpnEstimator only supports Designators and no composed expressions");
        auto spn_it = model->spns_.find(designator->table_name.text.assert_not_none());
        if (spn_it == model->spns_.end())
            throw data_model_exception("Could not find table for grouping.");

        auto &spn = spn_it->second.get();
        auto &attr_to_id = spn.get_attribute_to_id();
        if (auto attr_it = attr_to_id.find(designator->attr_name.text.assert_not_none()); attr_it != attr_to_id.end())
        {
            num_rows *= spn.estimate_number_distinct_values(attr_it->second);
        }
        else
        {
            num_rows *= spn.num_rows(); // if attribute is primary key, distinct values = num rows
        }
    }
    model->num_rows_ = num_rows;
    return model;
}

std::unique_ptr<DataModel>
SpnEstimator::estimate_join(const QueryGraph &, const DataModel &_left, const DataModel &_right,
                            const cnf::CNF &condition) const
{
    auto &left = as<const SpnDataModel>(_left);
    auto &right = as<const SpnDataModel>(_right);

    auto new_left = std::make_unique<SpnDataModel>(left);
    auto new_right = std::make_unique<SpnDataModel>(right);

    JoinTranslator jt;

    if (not condition.empty())
    {
        /* only consider single equi join */
        jt(*condition[0][0]);
        auto first_identifier = std::make_pair(jt.join_designator[0].first, jt.join_designator[0].second);
        auto second_identifier = std::make_pair(jt.join_designator[1].first, jt.join_designator[1].second);
        SpnJoin join = std::make_pair(first_identifier, second_identifier);

        /* if a data model is only responsible for one table (one spn) add the maximum frequency of the value
         * of the joined attribute */
        if (left.spns_.size() == 1)
        {
            new_left->max_frequencies_.emplace_back(max_frequency(left, join));
        }
        if (right.spns_.size() == 1)
        {
            new_right->max_frequencies_.emplace_back(max_frequency(right, join));
        }

        /* compute the estimated cardinality of the join via distinct count estimates.
         * See http://www.cidrdb.org/cidr2021/papers/cidr2021_paper01.pdf */
        std::size_t mf_left = std::accumulate(
            new_left->max_frequencies_.begin(), new_left->max_frequencies_.end(), 1, std::multiplies<>());

        std::size_t mf_right = std::accumulate(
            new_right->max_frequencies_.begin(), new_right->max_frequencies_.end(), 1, std::multiplies<>());

        std::size_t left_clause = new_left->num_rows_ / mf_left;
        std::size_t right_clause = new_right->num_rows_ / mf_right;

        std::size_t num_rows_join = std::min<std::size_t>(left_clause, right_clause) * mf_left * mf_right;

        new_left->num_rows_ = num_rows_join;
    }
    else
    {
        /* compute cartesian product since there is no join condition */
        if (left.spns_.size() == 1)
        {
            new_left->max_frequencies_.emplace_back(left.num_rows_);
        }
        if (right.spns_.size() == 1)
        {
            new_right->max_frequencies_.emplace_back(right.num_rows_);
        }

        new_left->num_rows_ = left.num_rows_ * right.num_rows_;
    }

    /* copy data from new_right to new_left to collect all information in one DataModel */
    new_left->spns_.insert(new_right->spns_.begin(), new_right->spns_.end());
    new_left->max_frequencies_.insert(
        new_left->max_frequencies_.end(), new_right->max_frequencies_.begin(), new_right->max_frequencies_.end());

    return new_left;
}

template <typename PlanTable>
std::unique_ptr<DataModel>
SpnEstimator::operator()(estimate_join_all_tag, PlanTable &&PT, const QueryGraph &, Subproblem to_join,
                         const cnf::CNF &condition) const
{
    M_insist(not to_join.empty());
    /* compute cartesian product */
    if (condition.empty())
    {
        auto model = std::make_unique<SpnDataModel>();
        model->num_rows_ = 1UL;
        for (auto it = to_join.begin(); it != to_join.end(); ++it)
            model->num_rows_ *= as<const SpnDataModel>(*PT[it.as_set()].model).num_rows_;
        return model;
    }

    /* get all attributes to join on */
    JoinTranslator jt;
    std::unordered_map<ThreadSafePooledString, ThreadSafePooledString> table_to_attribute;
    for (auto clause : condition)
    {
        jt(*clause[0]);
        table_to_attribute.emplace(jt.join_designator[0].first, jt.join_designator[0].second);
        table_to_attribute.emplace(jt.join_designator[1].first, jt.join_designator[1].second);
    }

    /* get first model to join */
    SpnDataModel final_model = as<const SpnDataModel>(*PT[to_join.begin().as_set()].model);
    ThreadSafePooledString first_table_name = final_model.spns_.begin()->first;

    /* if there is a join condition on this model, get the maximum frequency of the attribute */
    if (auto find_iter = table_to_attribute.find(first_table_name); find_iter != table_to_attribute.end())
    {
        final_model.max_frequencies_.emplace_back(max_frequency(final_model, find_iter->second));
    }
    /* else, maximum frequency is set as the number of rows (to compute the cartesian product) */
    else
    {
        final_model.max_frequencies_.emplace_back(final_model.spns_.begin()->second.get().num_rows());
    }

    /* iteratively add the next model to join via distinct count estimates */
    for (auto it = ++to_join.begin(); it != to_join.end(); it++)
    {
        SpnDataModel model = as<const SpnDataModel>(*PT[it.as_set()].model);
        ThreadSafePooledString table_name = model.spns_.begin()->first;

        if (auto find_iter = table_to_attribute.find(table_name); find_iter != table_to_attribute.end())
        {
            model.max_frequencies_.emplace_back(max_frequency(model, find_iter->second));
        }
        else
        {
            model.max_frequencies_.emplace_back(model.spns_.begin()->second.get().num_rows());
        }

        std::size_t mf_left = std::accumulate(
            final_model.max_frequencies_.begin(), final_model.max_frequencies_.end(), 1, std::multiplies<>());
        std::size_t mf_right = model.max_frequencies_[0];

        std::size_t left_clause = final_model.num_rows_ / mf_left;
        std::size_t right_clause = model.num_rows_ / mf_right;

        std::size_t num_rows_join = std::min<std::size_t>(left_clause, right_clause) * mf_left * mf_right;

        final_model.num_rows_ = num_rows_join;

        /* copy data from model to final_model to collect all information in one DataModel */
        final_model.spns_.emplace(*model.spns_.begin());
        final_model.max_frequencies_.emplace_back(model.max_frequencies_[0]);
    }

    return std::make_unique<SpnDataModel>(final_model);
}

std::size_t SpnEstimator::predict_cardinality(const DataModel &_data) const
{
    auto &data = as<const SpnDataModel>(_data);
    return data.num_rows_;
}

void SpnEstimator::print(std::ostream &) const {}
//=== SelectivityEstimator Constructors ===
SelectivityEstimator::SelectivityEstimator() = default;
SelectivityEstimator::SelectivityEstimator(ThreadSafePooledString name) { /* no-op */ }

//=== SelectivityEstimator operator() for join-all ===
template <typename PlanTable>
std::unique_ptr<DataModel>
SelectivityEstimator::operator()(estimate_join_all_tag, PlanTable &&PT, const QueryGraph &G, Subproblem to_join, const cnf::CNF &condition) const
{
    M_insist(!to_join.empty(), "SelectivityEstimator: to_join must not be empty");
    if (to_join.size() == 1)
    {
        auto &single_model = as<const SelectivityDataModel>(*PT[to_join.begin().as_set()].model);
        return std::make_unique<SelectivityDataModel>(single_model);
    }
    auto it = to_join.begin();
    auto &first_model = as<const SelectivityDataModel>(*PT[it.as_set()].model);
    auto result_model = std::make_unique<SelectivityDataModel>(first_model);
    ++it;
    for (; it != to_join.end(); ++it)
    {
        auto &right_model = as<const SelectivityDataModel>(*PT[it.as_set()].model);
        auto joined = estimate_join(G, *result_model, right_model, condition);
        result_model = std::unique_ptr<SelectivityDataModel>(static_cast<SelectivityDataModel *>(joined.release()));
    }
    return result_model;
}

//=== SelectivityEstimator predict_cardinality ===
std::size_t SelectivityEstimator::predict_cardinality(const DataModel &data) const
{
    return as<const SelectivityDataModel>(data).size;
}

//=== SelectivityEstimator print ===
void SelectivityEstimator::print(std::ostream &out) const
{
    out << "SelectivityEstimator - estimates cardinality using selectivity-based formulas";
}

std::unique_ptr<DataModel>
SelectivityEstimator::empty_model() const
{
    auto m = std::make_unique<SelectivityDataModel>();
    m->size = 0;
    return m;
}

std::unique_ptr<DataModel>
SelectivityEstimator::estimate_scan(const QueryGraph &G, Subproblem P) const
{
    M_insist(P.size() == 1, "Scan expects exactly one table");
    auto idx = *P.begin();
    auto &BT = as<const BaseTable>(*G.sources()[idx]);
    auto m = std::make_unique<SelectivityDataModel>();
    auto stats = BT.table().statistics();
    m->set_stats(*stats);
    m->size = stats->row_count;
    m->original_tables.insert(stats->table_name);
    return m;
}

std::unique_ptr<DataModel>
SelectivityEstimator::estimate_filter(const QueryGraph &G, const DataModel &data, const cnf::CNF &filter) const
{
    auto &dm = as<const SelectivityDataModel>(data);
    auto m = std::make_unique<SelectivityDataModel>(dm);
    if (filter.empty())
        return m;

    if (CardinalityStorage::Get().apply_stored_filter_cardinality(G, data, filter, *m)) {
        return m;
    }
    double sel = 1.0;
    auto stats = m->get_stats();
    auto filter_columns = filter.get_filter_columns();
    for (const auto &[table_name, column_set] : filter_columns)
    {
        for (const auto &column_name : column_set)
        {
            std::string col = table_name + "." + column_name;
            auto it = stats.selectivity.find(col);
            if (it != stats.selectivity.end())
                sel *= it->second;
        }
    }
    m->size = static_cast<std::size_t>(dm.size * sel);
    return m;
}

std::unique_ptr<DataModel>
SelectivityEstimator::estimate_join(const QueryGraph &G,
                                    const DataModel &left,
                                    const DataModel &right,
                                    const cnf::CNF &condition) const
{
    auto &lm = as<const SelectivityDataModel>(left);
    auto &rm = as<const SelectivityDataModel>(right);
    auto m = std::make_unique<SelectivityDataModel>();

    auto left_stats = lm.get_stats();
    auto right_stats = rm.get_stats();
    
    // Get all tables involved in left and right sides
    std::set<std::string> left_tables = lm.original_tables;
    std::set<std::string> right_tables = rm.original_tables;
    
    // Find join pairs that connect ANY table from left side to ANY table from right side
    std::vector<std::pair<std::string, std::string>> join_pairs;
    
    for (const auto &join : G.joins()) {
        auto join_condition = join->condition();
        if (join_condition.empty()) continue;
        
        auto join_columns = join_condition.get_join_columns();
        
        // Check if this join connects left side to right side
        bool connects_left_to_right = false;
        
        for (const auto &left_table : left_tables) {
            for (const auto &right_table : right_tables) {
                if (join_columns.count(left_table) && join_columns.count(right_table)) {
                    // Found a join that connects left_table to right_table
                    const auto &left_cols = join_columns.at(left_table);
                    const auto &right_cols = join_columns.at(right_table);
                    
                    if (!left_cols.empty() && !right_cols.empty()) {
                        std::string left_col = left_table + "." + *left_cols.begin();
                        std::string right_col = right_table + "." + *right_cols.begin();
                        join_pairs.emplace_back(left_col, right_col);
                        connects_left_to_right = true;
                    }
                }
            }
        }
        
        if (connects_left_to_right) break;
    }
    
    // Apply selectivity estimation
    double sel = 1.0;
    if (join_pairs.empty()) {
        m->size = lm.size * rm.size; // Cartesian product
    } else {
        for (const auto &[left_col, right_col] : join_pairs) {
            if (left_stats.distinct_counts.count(left_col) && 
                right_stats.distinct_counts.count(right_col)) {
                auto nd1 = left_stats.distinct_counts.at(left_col);
                auto nd2 = right_stats.distinct_counts.at(right_col);
                sel *= 1.0 / double(std::max(nd1, nd2));
            }
        }
        m->size = std::size_t(double(lm.size) * double(rm.size) * sel);
    }

    // Merge statistics and track original tables
    auto merged = left_stats.merge_for_join(right_stats);
    merged.row_count = m->size;
    
    // Combine original tables
    std::set<std::string> all_tables = left_tables;
    all_tables.insert(right_tables.begin(), right_tables.end());
    
    for (auto &kv : merged.distinct_counts) {
        merged.selectivity[kv.first] = double(kv.second) / double(merged.row_count);
    }
    
    m->set_stats(merged);
    m->original_tables = all_tables;
    return m;
}

std::unique_ptr<DataModel>
SelectivityEstimator::estimate_grouping(const QueryGraph &G, const DataModel &data, const std::vector<group_type> &groups) const
{
    auto &dm = as<const SelectivityDataModel>(data);
    auto m = std::make_unique<SelectivityDataModel>(dm);
    if (groups.empty())
    {
        m->size = 1;
        return m;
    }
    auto stats = m->get_stats();
    double prod = 1.0;
    
    for (auto &pr : groups)
    {
        std::string col = extract_column_name_from_expression(pr.first);
        
        // Use distinct_counts instead of histograms
        if (stats.distinct_counts.count(col)) {
            prod *= double(stats.distinct_counts.at(col));
        } else {
            // Fallback: assume column has moderate distinctness
            prod *= dm.size;
        }
    }
    
    m->size = static_cast<std::size_t>(prod);
    
    // Update statistics for the grouped result
    auto new_stats = m->get_stats();
    new_stats.row_count = m->size;
    
    // Recompute selectivities based on new row count
    for (auto &kv : new_stats.distinct_counts) {
        new_stats.selectivity[kv.first] = double(kv.second) / double(new_stats.row_count);
    }
    
    m->set_stats(new_stats);
    return m;
}

std::unique_ptr<DataModel>
SelectivityEstimator::estimate_limit(const QueryGraph &G, const DataModel &data, std::size_t limit, std::size_t offset) const
{
    auto &dm = as<const SelectivityDataModel>(data);
    auto m = std::make_unique<SelectivityDataModel>(dm);
    std::size_t rem = offset > dm.size ? 0 : dm.size - offset;
    m->size = std::min(rem, limit);
    return m;
}


std::string SelectivityEstimator::extract_column_name_from_expression(const ast::Expr &expr) const
{
    if (auto designator = cast<const ast::Designator>(&expr)) {
        if (designator->has_table_name()) {
            return std::string(*designator->table_name.text) + "." + *designator->attr_name.text;
        } else {
            return std::string(*designator->attr_name.text);
        }
    }
    return "";
}

#define LIST_CE(X)                                                                                   \
    X(CartesianProductEstimator, "CartesianProduct", "estimates cardinalities as Cartesian product") \
    X(InjectionCardinalityEstimator, "Injected", "estimates cardinalities based on a JSON file")     \
    X(SpnEstimator, "Spn", "estimates cardinalities based on Sum-Product Networks")                  \
    X(HistogramEstimator, "Histogram", "Used to test new estimation")                          \
    X(SelectivityEstimator, "Selectivitybased", "Used to test new estimation")                       \
    X(RangeCartesianProductEstimator, "RangeCartesianProduct", "range-aware cardinality estimator")

#define INSTANTIATE(TYPE, _1, _2)                                                                             \
    template std::unique_ptr<DataModel> TYPE::operator()(estimate_join_all_tag, PlanTableSmallOrDense &&PT,   \
                                                         const QueryGraph &G, Subproblem to_join,             \
                                                         const cnf::CNF &condition) const;                    \
    template std::unique_ptr<DataModel> TYPE::operator()(estimate_join_all_tag, PlanTableLargeAndSparse &&PT, \
                                                         const QueryGraph &G, Subproblem to_join,             \
                                                         const cnf::CNF &condition) const;
LIST_CE(INSTANTIATE)
#undef INSTANTIATE

__attribute__((constructor(202))) static void register_cardinality_estimators()
{
    Catalog &C = Catalog::Get();

#define REGISTER(TYPE, NAME, DESCRIPTION) \
    C.register_cardinality_estimator<TYPE>(C.pool(NAME), DESCRIPTION);
    LIST_CE(REGISTER)
#undef REGISTER

    C.arg_parser().add<bool>(
        /* group=       */ "Cardinality estimation",
        /* short=       */ nullptr,
        /* long=        */ "--show-cardinality-file-example",
        /* description= */ "show an example of an input JSON file for cardinality injection",
        [](bool)
        {
            std::cout << "\
Example for injected cardinalities file:\n\
{\n\
    database1: [\n\
            {\n\
                \"relations\": [\"A\", \"B\", ...],\n\
                \"size\": 150\n\
            },\n\
            {\n\
                \"relations\": [\"C\", \"A\", ...],\n\
                \"size\": 100\n\
            },\n\
    },\n\
    database2: [\n\
            {\n\
                \"relations\": [\"customers\"],\n\
                \"size\": 1000\n\
            },\n\
            {\n\
                \"relations\": [\"customers\", \"orders\", ...],\n\
                \"size\": 50\n\
            },\n\
    },\n\
}\n";
            exit(EXIT_SUCCESS);
        });
    C.arg_parser().add<const char *>(
        /* group=       */ "Cardinality estimation",
        /* short=       */ nullptr,
        /* long=        */ "--use-cardinality-file",
        /* description= */ "inject cardinalities from the given JSON file",
        [](const char *path)
        {
            options::injected_cardinalities_file = path;
        });
}
