#include <mutable/catalog/TableStatistics.hpp>
#include <mutable/catalog/Schema.hpp>
#include <mutable/IR/Tuple.hpp>
#include <backend/Interpreter.hpp>
#include <unordered_set>
#include <sstream>

namespace m
{

    void TableStatistics::compute(const Table &table)
    {
        row_count = table.store().num_rows();
        selectivity.clear();
        histograms.clear();

        const auto &schema = table.schema();
        Tuple tuple(schema);

        // Prepare a set for each column to track unique values
        std::vector<std::unordered_set<std::string>> unique_values(schema.num_entries()); // This is wrong, num_entries returns the size of the vector

        // For each row, use Interpreter to load the tuple
        for (std::size_t row = 0; row < row_count; ++row)
        {
            // Compile a loader for this row
            auto loader = Interpreter::compile_load(schema, table.store().memory().addr(), table.layout(), schema, row, 0);
            Tuple *args[] = {&tuple};
            loader(args);

            // For each column, collect unique values
            for (std::size_t col = 0; col < schema.num_entries(); ++col)
            {
                std::ostringstream oss;
                oss << tuple[col];
                unique_values[col].insert(oss.str());
            }
        }

        // Compute selectivity for each column
        for (std::size_t col = 0; col < schema.num_entries(); ++col)
        {
            double sel = row_count > 0 ? double(unique_values[col].size()) / double(row_count) : 1.0;
            // Use ThreadSafePooledString directly as the key
            selectivity[schema[col].id.name] = sel;
        }
    }

} // namespace m