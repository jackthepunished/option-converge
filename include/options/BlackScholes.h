#ifndef BLACK_SCHOLES_H
#define BLACK_SCHOLES_H

#include "PricingEngine.h"

namespace Options {

// Black-Scholes analytical pricing engine
class BlackScholes : public PricingEngine {
public:
    BlackScholes() = default;

    // Price an option using Black-Scholes formula
    PricingResult price(const OptionParams& params) override;

    // Calculate Greeks analytically
    Greeks calculateGreeks(const OptionParams& params) override;

    std::string getName() const override { return "Black-Scholes"; }

private:
    // Helper functions for Black-Scholes
    double d1(const OptionParams& params) const;
    double d2(const OptionParams& params) const;

    // Standard normal CDF
    double normalCDF(double x) const;

    // Standard normal PDF
    double normalPDF(double x) const;

    // Calculate individual option prices
    double callPrice(const OptionParams& params) const;
    double putPrice(const OptionParams& params) const;

    // Calculate individual Greeks
    double callDelta(const OptionParams& params) const;
    double putDelta(const OptionParams& params) const;
    double gamma(const OptionParams& params) const;
    double vega(const OptionParams& params) const;
    double callTheta(const OptionParams& params) const;
    double putTheta(const OptionParams& params) const;
    double callRho(const OptionParams& params) const;
    double putRho(const OptionParams& params) const;
};

} // namespace Options

#endif // BLACK_SCHOLES_H
