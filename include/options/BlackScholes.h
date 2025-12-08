#ifndef BLACK_SCHOLES_H
#define BLACK_SCHOLES_H

#include "PricingEngine.h"

namespace Options {

// Black-Scholes analytical pricing engine
class BlackScholes final : public PricingEngine {
public:
    BlackScholes() = default;
    ~BlackScholes() override = default;

    // Allow copying and moving
    BlackScholes(const BlackScholes&) = default;
    BlackScholes(BlackScholes&&) noexcept = default;
    BlackScholes& operator=(const BlackScholes&) = default;
    BlackScholes& operator=(BlackScholes&&) noexcept = default;

    // Price an option using Black-Scholes formula
    [[nodiscard]] PricingResult price(const OptionParams& params) override;

    // Calculate Greeks analytically
    [[nodiscard]] Greeks calculateGreeks(const OptionParams& params) override;

    [[nodiscard]] std::string getName() const noexcept override { return "Black-Scholes"; }

private:
    // Helper functions for Black-Scholes (inline for performance)
    [[nodiscard]] inline double d1(const OptionParams& params) const noexcept;
    [[nodiscard]] inline double d2(const OptionParams& params) const noexcept;

    // Standard normal CDF (inline for performance)
    [[nodiscard]] inline double normalCDF(double x) const noexcept;

    // Standard normal PDF (inline for performance)
    [[nodiscard]] inline double normalPDF(double x) const noexcept;

    // Calculate individual option prices (deprecated - kept for compatibility)
    [[nodiscard]] double callPrice(const OptionParams& params) const;
    [[nodiscard]] double putPrice(const OptionParams& params) const;

    // Calculate individual Greeks (deprecated - kept for compatibility)
    [[nodiscard]] double callDelta(const OptionParams& params) const;
    [[nodiscard]] double putDelta(const OptionParams& params) const;
    [[nodiscard]] double gamma(const OptionParams& params) const;
    [[nodiscard]] double vega(const OptionParams& params) const;
    [[nodiscard]] double callTheta(const OptionParams& params) const;
    [[nodiscard]] double putTheta(const OptionParams& params) const;
    [[nodiscard]] double callRho(const OptionParams& params) const;
    [[nodiscard]] double putRho(const OptionParams& params) const;
};

} // namespace Options

#endif // BLACK_SCHOLES_H
