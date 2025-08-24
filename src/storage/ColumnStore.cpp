#include "storage/ColumnStore.hpp"

#include "backend/StackMachine.hpp"
#include <mutable/catalog/Catalog.hpp>
#include <numeric>
#include <random>
#include <vector>
#include <cstring>
#include <backend/Interpreter.hpp>

using namespace m;

ColumnStore::ColumnStore(const Table &table)
    : Store(table)
{
    uint64_t max_attr_size = 0;

    /* Allocate memory for the attributes columns and the null bitmap column. */
    data_ = allocator_.allocate(ALLOCATION_SIZE * (table.num_attrs() + 1));

    /* Compute the capacity depending on the column with the largest attribute size. */
    for (auto attr = table.begin_all(); attr != table.end_all(); ++attr)
    {
        auto size = attr->type->size();
        row_size_ += size;
        max_attr_size = std::max(max_attr_size, size);
    }
    uint32_t num_attrs = table.num_attrs();
    row_size_ += num_attrs;
    uint64_t null_bitmap_size = /* pad_null_bitmap= */ 1 ? ((num_attrs + 7) / 8) * 8 : num_attrs;
    max_attr_size = std::max(max_attr_size, null_bitmap_size);

    capacity_ = (ALLOCATION_SIZE * 8) / max_attr_size;
}

ColumnStore::~ColumnStore() {}

// std::vector<Value> ColumnStore::sample_column(std::size_t attr_id, std::size_t sample_size, std::mt19937 &rng) const {
//     M_insist(attr_id < table().num_attrs(), "Attribute ID out of range");

//     const auto &attr_type = table().schema()[attr_id].type;

//     auto keep_test = table().schema();
//     auto test_two = keep_test.entries();
//     const Type *test_type = test_two[attr_id].type;
//     auto relevant_size = test_type->size() / 8;
//     auto aligned_relevant_size = test_type->alignment() / 8;

//     // Get the number of rows and the memory address of the column
//     std::size_t num_rows = this->num_rows();
//     uint8_t *column_data = reinterpret_cast<uint8_t *>(memory(attr_id));

//     std::vector<std::size_t> indices(num_rows);
//     std::iota(indices.begin(), indices.end(), 0);
//     std::shuffle(indices.begin(), indices.end(), rng);

//     std::vector<Value> sampled_values;
//     sampled_values.reserve(sample_size);

//     // Define a function pointer for extracting values
//     std::function<Value(const void *)> extract_value;

//     // Precompute the extraction logic based on the type
//     if (attr_type->is_integral()) {
//         if (relevant_size == 4) {
//             extract_value = [](const void *raw_data) -> Value {
//                 return Value(*static_cast<const int32_t *>(raw_data));
//             };
//         } else if (relevant_size == 8) {
//             extract_value = [](const void *raw_data) -> Value {
//                 return Value(*static_cast<const int64_t *>(raw_data));
//             };
//         } else {
//             throw std::runtime_error("Unsupported integral size");
//         }
//     } else if (attr_type->is_floating_point()) {
//         if (attr_type->is_float()) {
//             extract_value = [](const void *raw_data) -> Value {
//                 return Value(*static_cast<const float *>(raw_data));
//             };
//         } else if (attr_type->is_double()) {
//             extract_value = [](const void *raw_data) -> Value {
//                 return Value(*static_cast<const double *>(raw_data));
//             };
//         } else {
//             throw std::runtime_error("Unsupported floating-point size");
//         }
//     } else if (attr_type->is_boolean()) {
//         extract_value = [](const void *raw_data) -> Value {
//             return Value(*static_cast<const bool *>(raw_data));
//         };
//     } else {
//         throw std::runtime_error("Unsupported attribute type!");
//     }

//     // Use the precomputed extraction logic in the loop
//     for (std::size_t i = 0; i < sample_size; ++i) {
//         const void *raw_data = column_data + indices[i] * aligned_relevant_size;
//         sampled_values.push_back(extract_value(raw_data));
//         std::cout << "Value at index " << indices[i] << " is " << extract_value(raw_data) << std::endl;
//     }

//     return sampled_values;
// }

std::vector<Tuple> ColumnStore::sample_column(std::size_t attr_id, std::size_t sample_size, std::mt19937 &rng) const
{
    // Get the schema of the table
    const Schema &schema = table().schema();
    data_.dump(std::cout);

        // Get the number of rows and the memory address of the store
        std::size_t num_rows = this->num_rows();

    // Shuffle row indices for sampling
    std::vector<std::size_t> indices(num_rows);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    // Create a stack (or vector) to store sampled tuples
    std::vector<Tuple> sampled_tuples;
    sampled_tuples.reserve(sample_size);

    for (std::size_t row_id = 0; row_id < num_rows; ++row_id)
    {
        std::cout << "Row " << row_id << ": ";
        for (std::size_t col = 0; col < schema.num_entries(); ++col)
        {
            const auto &entry = schema[col];
            const Type *type = entry.type;

            uint8_t *store_data = reinterpret_cast<uint8_t *>(memory(col));
            std::size_t offset = row_id * (type->size() / 8);

            std::cout << "Column " << col << ": ";
            for (std::size_t i = 0; i < type->size() / 8; ++i)
            {
                std::cout << std::hex << static_cast<int>(store_data[offset + i]) << " ";
            }
            std::cout << std::dec << std::endl;
        }
    }

    std::vector<std::pair<void *, int32_t>> non_zero_values;

    // Start at the memory address of the first column
    uint8_t *start_address = reinterpret_cast<uint8_t *>(memory(0));
    uint8_t *end_address = start_address + capacity_ * row_size_ / 8; // Total memory size in bytes

    // Iterate over the memory space in 4-byte chunks
    for (uint8_t *current_address = start_address; current_address < end_address; current_address += 4)
    {
        int32_t value = *reinterpret_cast<const int32_t *>(current_address);

        // Store the address and value if the value is non-zero
        if (value != 0)
        {
            non_zero_values.emplace_back(static_cast<void *>(current_address), value);
        }
    }

    std::cout << "Non-zero values in memory:" << std::endl;

    void *last_address = nullptr; // To store the last address for calculating the delta

    for (const auto &[address, value] : non_zero_values) {
        std::cout << "Address: " << address << ", Value: " << value;

        // Calculate and print the delta to the last address
        if (last_address != nullptr) {
            auto delta = reinterpret_cast<uint8_t *>(address) - reinterpret_cast<uint8_t *>(last_address);
            std::cout << ", Delta to last: " << delta << " bytes";
        }

        std::cout << std::endl;

        // Update the last address
        last_address = address;
    }

    for (std::size_t col = 0; col < table().num_attrs(); ++col) {
        auto address = reinterpret_cast<uint8_t *>(memory(col));
        std::cout << "Column " << col << ": Base address = " << static_cast<void *>(address) << std::endl;
    }

    for (std::size_t col = 0; col < table().num_attrs(); ++col) {
        auto address = reinterpret_cast<uintptr_t>(memory(col));
        std::cout << "Column " << col << ": Address = " << address
                << ", Alignment = " << (address % 64) << " bytes" << std::endl;
    }

    // Sample tuples
    for (std::size_t i = 0; i < sample_size && i < num_rows; ++i)
    {
        std::size_t row_id = indices[i];

        // Create a new tuple based on the schema
        Tuple tuple(schema);

        // Fill the tuple with values from the store
        for (std::size_t col = 0; col < schema.num_entries(); ++col)
        {
            const auto &entry = schema[col];
            const Type *type = entry.type;

            uint8_t *store_data = reinterpret_cast<uint8_t *>(memory(col));
            std::size_t offset = row_id * (type->size() / 8);
            std::size_t offset_aligned = row_id * (type->alignment() / 8);

            if (type->is_integral())
            {
                if (type->size() / 8 == 4)
                {
                    int32_t value = *reinterpret_cast<const int32_t *>(store_data + offset);
                    tuple.set(col, value);
                }
                else if (type->size() / 8 == 8)
                {
                    int64_t value = *reinterpret_cast<const int64_t *>(store_data + offset);
                    tuple.set(col, value);
                }
            }
            else if (type->is_floating_point())
            {
                if (type->is_float())
                {
                    float value = *reinterpret_cast<const float *>(store_data + offset);
                    tuple.set(col, value);
                }
                else if (type->is_double())
                {
                    double value = *reinterpret_cast<const double *>(store_data + offset);
                    tuple.set(col, value);
                }
            }
            else if (type->is_character_sequence())
            {
                const char *value = reinterpret_cast<const char *>(store_data + offset);
                tuple.set(col, value);
            }
            else if (type->is_boolean())
            {
                bool value = *reinterpret_cast<const bool *>(store_data + offset);
                tuple.set(col, value);
            }
            else
            {
                throw std::runtime_error("Unsupported attribute type!");
            }
        }

        // Add the tuple to the stack
        tuple.print(std::cout, schema);
        std::cout << std::endl;
        sampled_tuples.push_back(std::move(tuple));
    }

    return sampled_tuples;
}

M_LCOV_EXCL_START
void ColumnStore::dump(std::ostream &out) const
{
    out << "ColumnStore for table \"" << table().name() << "\": " << num_rows_ << '/' << capacity_
        << " rows, " << row_size_ << " bits per row" << std::endl;
}
M_LCOV_EXCL_STOP

__attribute__((constructor(202))) static void register_store()
{
    Catalog &C = Catalog::Get();
    C.register_store<ColumnStore>(C.pool("ColumnStore"), "stores attributes in column-major order");
}
