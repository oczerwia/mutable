This is an excellent way to build a simple, yet principled, optimizer for benchmarking. Here are all the formulas you need for a selectivity-based estimator that operates in a single pass, using only basic statistics like table cardinality ($|T|$) and number of distinct values ($\text{n\_distinct}$).

All formulas rely on two key assumptions:
1.  **Uniformity:** Values are distributed uniformly across the data, unless explicitly captured by more detailed statistics.
2.  **Independence:** Different predicates and columns are independent of each other, allowing their selectivities to be multiplied.

---

### 1. Filter Predicate Operator

This operator estimates the number of rows remaining after a `WHERE` clause is applied. It takes the cardinality of the input relation as an argument.

$$\text{Card(Filter)} = \text{InputCard} \times \text{TotalSelectivity}$$

The `TotalSelectivity` is the product of the selectivities of each individual predicate:

$$\text{TotalSelectivity} = \text{selectivity}(P_1) \times \text{selectivity}(P_2) \times \dots$$

**Formulas for Individual Predicate Selectivity:**

* **Equality Predicate ($T.C = \text{value}$):**
    Assumes the value exists and all values are uniformly distributed.
    $$\text{selectivity}(T.C = \text{value}) = \frac{1}{\text{n\_distinct}(T.C)}$$

* **Inequality Predicate ($T.C < \text{value}$, $T.C > \text{value}$, etc.):**
    A simple and common heuristic for a prototype is to assume a fixed fraction of the data.
    $$\text{selectivity}(T.C < \text{value}) = \text{selectivity}(T.C > \text{value}) = \frac{1}{3}$$
    For `BETWEEN`, you could use a formula like $1/2$ or $1/3$, as it's a range.

* **`LIKE` Predicate ($T.C LIKE \text{'prefix\%'}$):**
    A common heuristic is to assume a fixed fraction of the data. For a prefix of length $L$, you could use a simple function like:
    $$\text{selectivity}(T.C \text{ LIKE 'prefix\%'}) = \frac{1}{100}$$
    (Or another small, fixed constant).

### 2. Group By Operator

This operator estimates the number of groups (i.e., rows in the final result set).

$$\text{Card(Group By)} = \text{min}(\text{InputCard}, \text{TotalNumGroups})$$

The `TotalNumGroups` is the product of the number of distinct values for each grouping key, based on the independence assumption.

$$\text{TotalNumGroups} = \text{n\_distinct}(G_1) \times \text{n\_distinct}(G_2) \times \dots$$

This formula provides a safe upper bound by taking the minimum of the calculated number of groups and the total number of rows in the input to the `GROUP BY` operator.

### 3. Join Operator (Multi-Join)

This operator estimates the cardinality of a multi-table join in a single pass. This is the formula you mentioned previously, restated here for completeness.

$$\text{Card(Join)} = (\text{InputCard}_1 \times \dots \times \text{InputCard}_n) \times \text{TotalJoinSelectivity}$$

**Important:** The `InputCard` for each table should be the cardinality *after* any `WHERE` clause predicates have been applied to it.

The `TotalJoinSelectivity` is the product of the selectivities of each individual join predicate.
$$\text{TotalJoinSelectivity} = \text{selectivity}(P_1) \times \text{selectivity}(P_2) \times \dots$$
where for each predicate:
$$\text{selectivity}(T_i.x = T_j.y) = \frac{1}{\text{max}(\text{n\_distinct}(T_i.x), \text{n\_distinct}(T_j.y))}$$