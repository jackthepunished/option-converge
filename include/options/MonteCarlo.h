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
class MonteCarlo final : public PricingEngine {
public:
    // Constructor
    explicit MonteCarlo(size_t numPaths = 10000,
                       DiscretizationScheme scheme = DiscretizationScheme::EULER,
                       VarianceReduction varRed = VarianceReduction::ANTITHETIC,
                       unsigned int seed = 42);
    ~MonteCarlo() override = default;

    // Allow copying and moving
    MonteCarlo(const MonteCarlo&) = default;
    MonteCarlo(MonteCarlo&&) noexcept = default;
    MonteCarlo& operator=(const MonteCarlo&) = default;
    MonteCarlo& operator=(MonteCarlo&&) noexcept = default;

    // Price an option using Monte Carlo simulation
    [[nodiscard]] PricingResult price(const OptionParams& params) override;

    // Calculate Greeks using pathwise derivatives or finite differences
    [[nodiscard]] Greeks calculateGreeks(const OptionParams& params) override;

    [[nodiscard]] std::string getName() const override;

    // Setters
    void setNumPaths(size_t numPaths) noexcept { numPaths_ = numPaths; }
    void setScheme(DiscretizationScheme scheme) noexcept { scheme_ = scheme; }
    void setVarianceReduction(VarianceReduction varRed) noexcept { varRed_ = varRed; }
    void setSeed(unsigned int seed);

    // Getters
    [[nodiscard]] size_t getNumPaths() const noexcept { return numPaths_; }
    [[nodiscard]] DiscretizationScheme getScheme() const noexcept { return scheme_; }
    [[nodiscard]] VarianceReduction getVarianceReduction() const noexcept { return varRed_; }

private:
    size_t numPaths_;
    DiscretizationScheme scheme_;
    VarianceReduction varRed_;
    std::mt19937 rng_;

    // Simulate one path using Euler scheme
    [[nodiscard]] double simulatePathEuler(const OptionParams& params);

    // Simulate one path using Milstein scheme
    [[nodiscard]] double simulatePathMilstein(const OptionParams& params);

    // Calculate payoff at maturity
    [[nodiscard]] double payoff(double finalPrice, const OptionParams& params) const noexcept;

    // Standard normal random number
    [[nodiscard]] double randomNormal();

    // Price with basic Monte Carlo
    [[nodiscard]] double priceBasic(const OptionParams& params, double& stdError);

    // Price with antithetic variates
    [[nodiscard]] double priceAntithetic(const OptionParams& params, double& stdError);

    // Price with control variates (using Black-Scholes)
    [[nodiscard]] double priceControlVariate(const OptionParams& params, double& stdError);

    // Price with both variance reduction techniques
    [[nodiscard]] double priceBoth(const OptionParams& params, double& stdError);

    // Calculate standard error
    [[nodiscard]] double calculateStandardError(const std::vector<double>& payoffs) const;
};

} // namespace Options

#endif // MONTE_CARLO_H
