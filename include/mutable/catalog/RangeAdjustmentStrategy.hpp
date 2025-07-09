#include <utility>
#include <cmath>

struct RangeAdjustmentStrategy
{
    virtual ~RangeAdjustmentStrategy() = default;
    virtual std::pair<double, double> adjust(
        std::pair<double, double> estimated_range,
        double true_cardinality) const = 0;
};

// Example: Pull up lower bound by 10% of the error margin, round bounds
struct TightenBoundsStrategy : public RangeAdjustmentStrategy
{
    double factor; // e.g., 0.1 for 10%
    TightenBoundsStrategy(double factor = 0.3) : factor(factor) {}

    std::pair<double, double> adjust(
        std::pair<double, double> estimated_range,
        double true_cardinality) const override
    {
        double lower = estimated_range.first;
        double upper = estimated_range.second;
        // Case 1: Our true cardinality is within our range
        // We generously update our lower and upper bound
        if (true_cardinality >= lower && true_cardinality <= upper)
        {
            double lower_error_margin = true_cardinality - lower;
            lower = lower + factor * lower_error_margin;

            double upper_error_margin = upper - true_cardinality;
            upper = upper - factor * upper_error_margin;
        
        }
        // Case 2: Our true cardinality is below our range
        if (true_cardinality < lower && true_cardinality < upper)
                {
                    lower = lower - factor * lower;
                }
        // Case 3: Our true cardinality is above our range
        // This case will not be implemeted, as we implement a strict upper bound or use the cartesian product

        // Round as desired (here: up for lower, down for upper)
        return {std::ceil(lower), std::floor(upper)};
    }
};