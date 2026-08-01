#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "PricingEngine.h"

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
    unsigned int seed_;

    // Path work is partitioned into a fixed set of chunks, each with its own
    // deterministically seeded RNG, and the chunk loop is parallelised with
    // OpenMP where available. The partition - and therefore every random
    // draw - is independent of the thread count; estimates agree across
    // thread counts (and OpenMP-free builds) up to floating-point summation
    // order, around 1e-12 relative. Pricing is a pure function of the seed,
    // which is also what makes common-random-numbers Greeks work without
    // explicit RNG state management.

    // Calculate payoff at maturity
    [[nodiscard]] double payoff(double finalPrice, const OptionParams& params) const noexcept;

    // Price with basic Monte Carlo
    [[nodiscard]] double priceBasic(const OptionParams& params, double& stdError) const;

    // Price with antithetic variates
    [[nodiscard]] double priceAntithetic(const OptionParams& params, double& stdError) const;

    // Price with control variates (using the discounted terminal price)
    [[nodiscard]] double priceControlVariate(const OptionParams& params, double& stdError) const;

    // Price with both variance reduction techniques
    [[nodiscard]] double priceBoth(const OptionParams& params, double& stdError) const;
};

} // namespace Options

#endif // MONTE_CARLO_H
