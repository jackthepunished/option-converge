#ifndef IMPLIED_VOLATILITY_H
#define IMPLIED_VOLATILITY_H

#include "Option.h"

#include <cstddef>

namespace Options {

// Solves for the volatility that reproduces an observed market price under
// Black-Scholes. Newton-Raphson on the analytic vega converges quadratically
// near the money; where vega is too small to trust (deep in/out of the money,
// short maturities) the solver falls back to a bracketed Brent search, which
// is guaranteed to converge because the Black-Scholes price is strictly
// increasing in volatility.
class ImpliedVolatility {
public:
    struct Result {
        double impliedVol;    // The solved volatility (annualised, absolute)
        size_t iterations;    // Newton iterations used (excludes fallback)
        bool usedFallback;    // True if Brent finished the job
    };

    explicit ImpliedVolatility(double tolerance = 1e-8, size_t maxIterations = 100)
        : tolerance_(tolerance), maxIterations_(maxIterations) {}

    // Solve for the volatility at which Black-Scholes reproduces marketPrice.
    // The volatility field of params is ignored. Throws std::invalid_argument
    // if the option is not European or the price violates no-arbitrage bounds
    // (no volatility can reproduce it); throws std::runtime_error if the
    // solver fails to converge.
    [[nodiscard]] Result solve(double marketPrice, const OptionParams& params) const;

private:
    double tolerance_;
    size_t maxIterations_;

    // Black-Scholes price and vega (per unit volatility, not per 1%) at a
    // trial sigma. Free of PricingResult overhead so the solver loop is cheap.
    [[nodiscard]] static double priceAt(double sigma, const OptionParams& params) noexcept;
    [[nodiscard]] static double vegaAt(double sigma, const OptionParams& params) noexcept;

    [[nodiscard]] Result brent(double marketPrice, const OptionParams& params,
                               double lo, double hi, size_t newtonIterations) const;
};

} // namespace Options

#endif // IMPLIED_VOLATILITY_H
