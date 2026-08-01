#ifndef PRICING_ENGINE_H
#define PRICING_ENGINE_H

#include "Option.h"

namespace Options {

// Abstract base class for all pricing engines
class PricingEngine {
public:
    virtual ~PricingEngine() = default;

    // Pure virtual function for pricing
    [[nodiscard]] virtual PricingResult price(const OptionParams& params) = 0;

    // Calculate Greeks (default implementation uses finite differences)
    [[nodiscard]] virtual Greeks calculateGreeks(const OptionParams& params);

    // Get engine name
    [[nodiscard]] virtual std::string getName() const = 0;

protected:
    // Default constructor for derived classes
    PricingEngine() = default;

    // Copying through a base pointer or reference would slice, so the copy
    // and move operations are protected rather than public - and defaulted
    // rather than deleted, so that concrete engines, which are
    // self-contained value types, can default their own copy and move.
    // Deleting them here would silently delete the derived ones too.
    PricingEngine(const PricingEngine&) = default;
    PricingEngine(PricingEngine&&) = default;
    PricingEngine& operator=(const PricingEngine&) = default;
    PricingEngine& operator=(PricingEngine&&) = default;

    // Helper function for finite difference Greeks
    [[nodiscard]] double calculateDelta(const OptionParams& params, double h = 0.01);
    [[nodiscard]] double calculateGamma(const OptionParams& params, double h = 0.01);
    [[nodiscard]] double calculateVega(const OptionParams& params, double h = 0.01);
    [[nodiscard]] double calculateTheta(const OptionParams& params, double h = 1.0/365.0);
    [[nodiscard]] double calculateRho(const OptionParams& params, double h = 0.01);
};

} // namespace Options

#endif // PRICING_ENGINE_H
