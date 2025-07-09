// Used to define strategies of comparing ranges with each other for the plan enumerator

#include <utility>

struct RangeComparer
{
    // Simple interface to compare ranges - returns true if a is better than b
    virtual bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const = 0;
    virtual ~RangeComparer() = default;
};

// Concrete implementations
struct UpperBoundComparer : public RangeComparer
{
    bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const override
    {
        return a.second < b.second;
    }
};

struct LowerBoundComparer : public RangeComparer
{
    bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const override
    {
        return a.first < b.first;
    }
};


struct MeanUncertaintyComparer : public RangeComparer
{
    // weight: 0.0 = only mean, 1.0 = only uncertainty, in between = trade-off
    double uncertainty_weight;
    MeanUncertaintyComparer(double uncertainty_weight = 0.3)
        : uncertainty_weight(uncertainty_weight) {}

    bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const override
    {
        double mean_a = (a.first + a.second) / 2.0;
        double mean_b = (b.first + b.second) / 2.0;
        double unc_a = a.second - a.first;
        double unc_b = b.second - b.first;

        // Lower score is better
        double score_a = (1.0 - uncertainty_weight) * mean_a + uncertainty_weight * unc_a;
        double score_b = (1.0 - uncertainty_weight) * mean_b + uncertainty_weight * unc_b;

        return score_a < score_b;
    }
};