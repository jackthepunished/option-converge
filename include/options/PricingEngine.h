#ifndef PRICING_ENGINE_H
#define PRICING_ENGINE_H

#include "Option.h"

namespace Options {

// Abstract base class for all pricing engines
class PricingEngine {
public:
    virtual ~PricingEngine() = default;

    // Pure virtual function for pricing
    virtual PricingResult price(const OptionParams& params) = 0;

    // Calculate Greeks (default implementation uses finite differences)
    virtual Greeks calculateGreeks(const OptionParams& params);

    // Get engine name
    virtual std::string getName() const = 0;

protected:
    // Helper function for finite difference Greeks
    double calculateDelta(const OptionParams& params, double h = 0.01);
    double calculateGamma(const OptionParams& params, double h = 0.01);
    double calculateVega(const OptionParams& params, double h = 0.01);
    double calculateTheta(const OptionParams& params, double h = 1.0/365.0);
    double calculateRho(const OptionParams& params, double h = 0.01);
};

} // namespace Options

#endif // PRICING_ENGINE_H
