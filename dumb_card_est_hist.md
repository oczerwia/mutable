### Implementation Summary for a 1D Histogram-Based Cardinality Estimator

This report outlines the formulas and assumptions necessary to implement a cardinality estimator that relies exclusively on 1D equi-width histograms, distinct value counts (NDV), and total row counts. This model assumes independence between columns and uniform distribution of values within each histogram bucket.

#### 1. Histogram Structure and Assumptions

The core component is the 1D equi-width histogram. For a numeric column, the histogram is composed of a fixed number of buckets, where each bucket covers an equal-sized range of values.

* **Histogram Data Structure:** Each bucket $B_i$ contains the following information:
    * $low_i$: The minimum value in the bucket's range.
    * $high_i$: The maximum value in the bucket's range.
    * $count_i$: The number of rows in the table whose value falls within the bucket's range.
* **Fundamental Assumptions:**
    * **Uniformity Assumption:** Values (and distinct values) are uniformly distributed within each bucket. If a value is sought within a bucket, its estimated frequency is the total count for that bucket divided by the number of distinct values in the bucket's range.
    * **Independence Assumption:** For multiple predicates on different columns or for join conditions, the selectivities are independent and can be multiplied.

#### 2. Cardinality Estimation for Operators

Let's assume a table $T$ with `Total_Rows` and a column $C$. The histogram for $C$ is composed of $N$ buckets, $B_1, ..., B_N$. The total number of distinct values for column $C$ is $NDV_C$.

##### a) Filter Predicates

**Equality Predicate:** $WHERE \ C = V$
1.  Locate the bucket $B_i$ that contains the value $V$.
2.  Under the uniformity assumption, the estimated frequency of any value in bucket $B_i$ is its row count divided by the number of distinct values. Since we don't store distinct values per bucket, we assume the distinct values are uniformly distributed, giving us an estimated frequency. The most common heuristic is to use the bucket's width as the number of distinct values.
3.  The estimated cardinality is the number of rows in the bucket divided by its width.

* **Formula:**
    $Cardinality = \frac{count_i}{high_i - low_i}$

**Range/Inequality Predicates:** $WHERE \ C > V$ or $WHERE \ V_1 < C < V_2$
1.  Identify all buckets that fall entirely within the specified range.
2.  For any buckets that are only partially covered by the range, assume uniformity within the bucket to calculate the fraction of rows covered.
3.  Sum the counts of the fully covered buckets and the fractional counts of the partially covered buckets.

* **Formula for $WHERE \ V_1 < C < V_2$:**
    $Cardinality = \sum_{B_i \text{ fully in range}} count_i + \sum_{B_j \text{ partially in range}} count_j \cdot \frac{\text{covered_range_in_} B_j}{\text{bucket_width_of_} B_j}$

**Multiple Filter Predicates:** $WHERE \ C_1 = V_1 \ AND \ C_2 < V_2$
1.  Calculate the selectivity for each predicate individually, using the methods above.
2.  Multiply the individual selectivities due to the independence assumption.
3.  Multiply the combined selectivity by the total number of rows.

* **Formula:**
    $Selectivity(Total) = Selectivity(C_1 = V_1) \cdot Selectivity(C_2 < V_2)$
    $Cardinality = Total\_Rows \cdot Selectivity(Total)$

##### b) Equi-Join

For an equi-join between two tables $R$ and $S$ on columns $A$ and $B$, respectively ($R.A = S.B$), we must align their histograms.

1.  Identify the overlapping buckets between the histograms of $R.A$ and $S.B$.
2.  For each overlapping bucket, multiply the row counts from both tables for that bucket.
3.  Sum these products over all overlapping buckets.

* **Formula:**
    $Cardinality(R \Join S) = \sum_{i=1}^{N_{bins}} \frac{count_R(B_i) \cdot count_S(B_i)}{\text{bucket_width_} i}$

##### c) Multi-Join

Multi-joins (e.g., $R \Join S \Join T$) are estimated by sequentially applying the two-way join formula.

1.  First, estimate the cardinality of the join of the first two tables, say $R$ and $S$, using the equi-join formula above.
2.  The resulting cardinality is treated as the size of a temporary intermediate relation.
3.  Next, perform the join between this intermediate relation and the third table, $T$. The cardinality is the cardinality of the intermediate relation multiplied by the selectivity of the second join.
4.  The independence assumption is key, as it allows us to multiply join selectivities.

##### d) GROUP BY

The cardinality of a `GROUP BY` operation is the number of distinct groups. This is equivalent to estimating the number of distinct values (NDV) in the grouping column.

1.  If no filter is applied, the cardinality is simply the total distinct value count of the column.

* **Formula:**
    $Cardinality(GROUP \ BY \ C) = NDV_C$

2.  If the `GROUP BY` follows a filter, the estimated cardinality is the number of distinct values within the filtered subset. Without distinct values per bucket, a simple heuristic is to take the minimum of the filtered cardinality and the total distinct value count of the column. This prevents overestimating the number of groups beyond what is possible.

* **Formula:**
    $Cardinality(GROUP \ BY \ C \ after \ Filter) = min(NDV_C, \ Cardinality_{after\_filter})$