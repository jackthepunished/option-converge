#ifndef PRICING_ENGINE_H
#define PRICING_ENGINE_H

#include "Option.h"

namespace Options {

// Abstract base class for all pricing engines
class PricingEngine {
public:
    virtual ~PricingEngine() = default;

    // Prevent copying and moving of base class
    PricingEngine(const PricingEngine&) = delete;
    PricingEngine(PricingEngine&&) = delete;
    PricingEngine& operator=(const PricingEngine&) = delete;
    PricingEngine& operator=(PricingEngine&&) = delete;

    // Pure virtual function for pricing
    [[nodiscard]] virtual PricingResult price(const OptionParams& params) = 0;

    // Calculate Greeks (default implementation uses finite differences)
    [[nodiscard]] virtual Greeks calculateGreeks(const OptionParams& params);

    // Get engine name
    [[nodiscard]] virtual std::string getName() const = 0;

protected:
    // Default constructor for derived classes
    PricingEngine() = default;

    // Helper function for finite difference Greeks
    [[nodiscard]] double calculateDelta(const OptionParams& params, double h = 0.01);
    [[nodiscard]] double calculateGamma(const OptionParams& params, double h = 0.01);
    [[nodiscard]] double calculateVega(const OptionParams& params, double h = 0.01);
    [[nodiscard]] double calculateTheta(const OptionParams& params, double h = 1.0/365.0);
    [[nodiscard]] double calculateRho(const OptionParams& params, double h = 0.01);
};

} // namespace Options

#endif // PRICING_ENGINE_H
