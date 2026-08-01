#include "options/ImpliedVolatility.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Options {

namespace {

constexpr double kTwoPi = 6.283185307179586476925287;
constexpr double kSigmaLo = 1e-6;   // Below this the price is numerically flat
constexpr double kSigmaHi = 5.0;    // 500% vol: above any quoted market

double normalCDF(double x) noexcept {
    constexpr double kInvSqrt2 = 0.70710678118654752440;
    return 0.5 * std::erfc(-x * kInvSqrt2);
}

double normalPDF(double x) noexcept {
    constexpr double kInvSqrt2Pi = 0.3989422804014327;
    return kInvSqrt2Pi * std::exp(-0.5 * x * x);
}

} // namespace

double ImpliedVolatility::priceAt(double sigma, const OptionParams& params) noexcept {
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double T = params.timeToMaturity;
    const double sqrtT = std::sqrt(T);

    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;
    const double discS = S * std::exp(-q * T);
    const double discK = K * std::exp(-r * T);

    return params.isCall() ? discS * normalCDF(d1) - discK * normalCDF(d2)
                           : discK * normalCDF(-d2) - discS * normalCDF(-d1);
}

double ImpliedVolatility::vegaAt(double sigma, const OptionParams& params) noexcept {
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double T = params.timeToMaturity;
    const double sqrtT = std::sqrt(T);

    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    return S * std::exp(-q * T) * normalPDF(d1) * sqrtT;
}

ImpliedVolatility::Result ImpliedVolatility::solve(double marketPrice,
                                                   const OptionParams& params) const {
    if (!params.isEuropean()) {
        throw std::invalid_argument(
            "Implied volatility is Black-Scholes based and requires a European option.");
    }

    // No volatility can push the price outside the no-arbitrage band
    // (discounted intrinsic, discounted spot-or-strike); reject rather than
    // letting the solver chase a root that does not exist.
    const double discS = params.spotPrice * std::exp(-params.dividendYield * params.timeToMaturity);
    const double discK = params.strikePrice * std::exp(-params.riskFreeRate * params.timeToMaturity);
    const double lower = std::max(0.0, params.isCall() ? discS - discK : discK - discS);
    const double upper = params.isCall() ? discS : discK;
    if (marketPrice <= lower || marketPrice >= upper) {
        throw std::invalid_argument(
            "Market price violates no-arbitrage bounds; no implied volatility exists.");
    }

    // Brenner-Subrahmanyam starting point: exact for at-the-money-forward
    // options, and a good enough neighbourhood everywhere Newton is usable.
    double sigma = std::clamp(
        std::sqrt(kTwoPi / params.timeToMaturity) * marketPrice / params.spotPrice,
        0.05, 2.0);

    for (size_t i = 0; i < maxIterations_; ++i) {
        const double priceError = priceAt(sigma, params) - marketPrice;
        if (std::fabs(priceError) < tolerance_) {
            return {sigma, i, false};
        }
        const double vega = vegaAt(sigma, params);
        // A flat price surface makes the Newton step huge and unreliable —
        // exactly the low-vega regime the Brent fallback exists for.
        if (vega < 1e-10) {
            return brent(marketPrice, params, kSigmaLo, kSigmaHi, i);
        }
        const double next = sigma - priceError / vega;
        if (next <= kSigmaLo || next >= kSigmaHi) {
            return brent(marketPrice, params, kSigmaLo, kSigmaHi, i);
        }
        sigma = next;
    }
    return brent(marketPrice, params, kSigmaLo, kSigmaHi, maxIterations_);
}

// Classic Brent root-finding on f(sigma) = price(sigma) - marketPrice.
// Monotonicity of the Black-Scholes price in volatility guarantees the
// bracket [lo, hi] holds exactly one root once the bounds check has passed.
ImpliedVolatility::Result ImpliedVolatility::brent(double marketPrice,
                                                   const OptionParams& params,
                                                   double lo, double hi,
                                                   size_t newtonIterations) const {
    double a = lo;
    double b = hi;
    double fa = priceAt(a, params) - marketPrice;
    double fb = priceAt(b, params) - marketPrice;
    if (fa * fb > 0.0) {
        throw std::runtime_error(
            "Implied volatility root is not bracketed; price is at the edge of "
            "the representable range.");
    }

    double c = a;
    double fc = fa;
    double d = b - a;
    double e = d;

    for (size_t iter = 0; iter < 200; ++iter) {
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        const double tol = 2.0 * 2.220446049250313e-16 * std::fabs(b) + 0.5 * tolerance_;
        const double m = 0.5 * (c - b);
        if (std::fabs(m) <= tol || fb == 0.0) {
            return {b, newtonIterations, true};
        }
        if (std::fabs(e) < tol || std::fabs(fa) <= std::fabs(fb)) {
            d = m; e = m;  // Bisection
        } else {
            double p, q_;
            double s = fb / fa;
            if (a == c) {
                p = 2.0 * m * s;          // Secant
                q_ = 1.0 - s;
            } else {                      // Inverse quadratic interpolation
                const double qq = fa / fc;
                const double rr = fb / fc;
                p = s * (2.0 * m * qq * (qq - rr) - (b - a) * (rr - 1.0));
                q_ = (qq - 1.0) * (rr - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q_ = -q_;
            p = std::fabs(p);
            if (2.0 * p < std::min(3.0 * m * q_ - std::fabs(tol * q_), std::fabs(e * q_))) {
                e = d; d = p / q_;        // Accept interpolation
            } else {
                d = m; e = m;             // Fall back to bisection
            }
        }
        a = b;
        fa = fb;
        b += (std::fabs(d) > tol) ? d : (m > 0.0 ? tol : -tol);
        fb = priceAt(b, params) - marketPrice;
        if ((fb > 0.0) == (fc > 0.0)) {
            c = a; fc = fa;
            d = b - a; e = d;
        }
    }
    throw std::runtime_error("Implied volatility solver failed to converge.");
}

} // namespace Options
