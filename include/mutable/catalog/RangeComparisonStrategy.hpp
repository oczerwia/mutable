// Used to define strategies of comparing ranges with each other for the plan enumerator

#include <utility>
#include <mutable/Options.hpp>
#include <memory>

struct RangeComparer
{
    // Simple interface to compare ranges - returns true if a is better than b
    virtual bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const = 0;
    virtual ~RangeComparer() = default;
    virtual double collapse(const std::pair<double, double> &a) = 0;
};

// Concrete implementations
struct UpperBoundComparer : public RangeComparer
{
    bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const override
    {
        return a.second < b.second;
    }
    double collapse(const std::pair<double, double> &a){
        return a.second;
    }

};

struct LowerBoundComparer : public RangeComparer
{
    bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const override
    {
        return a.first < b.first;
    }
    double collapse(const std::pair<double, double> &a){
        return a.first;
    }
};


struct MeanComparer : public RangeComparer
{
    bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const override
    {
        double mean_a = (a.first + a.second) / 2.0;
        double mean_b = (b.first + b.second) / 2.0;
        return mean_a < mean_b;
    }
    double collapse(const std::pair<double, double> &a) override
    {
        return (a.first + a.second) / 2.0;
    }
};

struct MeanUncertaintyComparer : public RangeComparer
{
    bool compare(const std::pair<double, double> &a, const std::pair<double, double> &b) const override
    {
        double uncertainty_weight = m::Options::Get().uncertainty_impact;
        double mean_a = (a.first + a.second) / 2.0;
        double mean_b = (b.first + b.second) / 2.0;
        double unc_a = a.second - a.first;
        double unc_b = b.second - b.first;

        double score_a = (1.0 - uncertainty_weight) * mean_a + uncertainty_weight * unc_a;
        double score_b = (1.0 - uncertainty_weight) * mean_b + uncertainty_weight * unc_b;

        return score_a < score_b;
    }
    double collapse(const std::pair<double, double> &a){
        double uncertainty_weight = m::Options::Get().uncertainty_impact;
        double mean = (a.first + a.second) / 2.0;
        double unc = a.second - a.first;
        double score = (1.0 - uncertainty_weight) * mean + uncertainty_weight * unc;
        return score;
    }
};

inline std::unique_ptr<RangeComparer> GetRangeComparer_() {
    const char* strategy = m::Options::Get().collapse_function;
    if (!strategy || std::string(strategy) == "UpperBound") {
        return std::make_unique<UpperBoundComparer>();
    } else if (std::string(strategy) == "LowerBound") {
        return std::make_unique<LowerBoundComparer>();
    } else if (std::string(strategy) == "Uncertainty") {
        return std::make_unique<MeanUncertaintyComparer>();
    } else if (std::string(strategy) == "Mean") {
        return std::make_unique<MeanComparer>();
    } else {
        return std::make_unique<UpperBoundComparer>();
    }
}