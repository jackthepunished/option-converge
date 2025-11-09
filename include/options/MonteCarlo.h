#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "PricingEngine.h"
#include <random>
#include <vector>

namespace Options {

// Enum for discretization schemes
enum class DiscretizationScheme {
    EULER,
    MILSTEIN
};

// Enum for variance reduction techniques
enum class VarianceReduction {
    NONE,
    ANTITHETIC,
    CONTROL_VARIATE,
    BOTH
};

// Monte Carlo simulation pricing engine
class MonteCarlo : public PricingEngine {
public:
    // Constructor
    explicit MonteCarlo(size_t numPaths = 10000,
                       DiscretizationScheme scheme = DiscretizationScheme::EULER,
                       VarianceReduction varRed = VarianceReduction::ANTITHETIC,
                       unsigned int seed = 42);

    // Price an option using Monte Carlo simulation
    PricingResult price(const OptionParams& params) override;

    // Calculate Greeks using pathwise derivatives or finite differences
    Greeks calculateGreeks(const OptionParams& params) override;

    std::string getName() const override;

    // Setters
    void setNumPaths(size_t numPaths) { numPaths_ = numPaths; }
    void setScheme(DiscretizationScheme scheme) { scheme_ = scheme; }
    void setVarianceReduction(VarianceReduction varRed) { varRed_ = varRed; }
    void setSeed(unsigned int seed);

    // Getters
    size_t getNumPaths() const { return numPaths_; }
    DiscretizationScheme getScheme() const { return scheme_; }
    VarianceReduction getVarianceReduction() const { return varRed_; }

private:
    size_t numPaths_;
    DiscretizationScheme scheme_;
    VarianceReduction varRed_;
    std::mt19937 rng_;

    // Simulate one path using Euler scheme
    double simulatePathEuler(const OptionParams& params);

    // Simulate one path using Milstein scheme
    double simulatePathMilstein(const OptionParams& params);

    // Calculate payoff at maturity
    double payoff(double finalPrice, const OptionParams& params) const;

    // Standard normal random number
    double randomNormal();

    // Price with basic Monte Carlo
    double priceBasic(const OptionParams& params, double& stdError);

    // Price with antithetic variates
    double priceAntithetic(const OptionParams& params, double& stdError);

    // Price with control variates (using Black-Scholes)
    double priceControlVariate(const OptionParams& params, double& stdError);

    // Price with both variance reduction techniques
    double priceBoth(const OptionParams& params, double& stdError);

    // Calculate standard error
    double calculateStandardError(const std::vector<double>& payoffs) const;
};

} // namespace Options

#endif // MONTE_CARLO_H
