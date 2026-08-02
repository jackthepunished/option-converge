#include "options/LookbackOption.h"

#include <cmath>
#include <random>
#include <stdexcept>

namespace Options {

namespace {

double normalCDF(double x) noexcept {
    constexpr double kInvSqrt2 = 0.70710678118654752440;
    return 0.5 * std::erfc(-x * kInvSqrt2);
}

} // namespace

// Goldman-Sosin-Gatto (floating) and Conze-Viswanathan (fixed) prices at
// inception, where the running extreme equals the spot. b = r - q is the
// cost of carry; the (S/x)^(-2b/sigma^2) reflection terms simplify where
// the extreme is the spot itself.
double AnalyticLookback::price(const LookbackParams& params) {
    const OptionParams& o = params.option;
    const double S = o.spotPrice;
    const double K = o.strikePrice;
    const double r = o.riskFreeRate;
    const double q = o.dividendYield;
    const double sigma = o.volatility;
    const double T = o.timeToMaturity;
    const double b = r - q;

    if (std::fabs(b) < 1e-8) {
        throw std::invalid_argument(
            "Lookback closed forms carry a sigma^2/(2(r-q)) factor; the r == q limit "
            "is not implemented. Use the Monte Carlo pricer for that regime.");
    }

    const double sqrtT = std::sqrt(T);
    const double sigT = sigma * sqrtT;
    const double coef = sigma * sigma / (2.0 * b);
    const double reflT = (2.0 * b / sigma) * sqrtT;  // (2b/sigma) * sqrt(T)
    const double discR = std::exp(-r * T);
    const double carry = std::exp((b - r) * T);      // = e^{-qT}
    const double ebT = std::exp(b * T);

    // At inception ln(S/extreme) = 0, so the at-the-extreme d-value is:
    const double a1 = (b + 0.5 * sigma * sigma) * T / sigT;
    const double a2 = a1 - sigT;

    if (params.lookbackType == LookbackType::FLOATING_STRIKE) {
        if (o.isCall()) {
            // Pays S_T - min(S): always in the money at expiry.
            return S * carry * normalCDF(a1) - S * discR * normalCDF(a2) +
                   S * discR * coef *
                       (normalCDF(-a1 + reflT) - ebT * normalCDF(-a1));
        }
        // Pays max(S) - S_T.
        return S * discR * normalCDF(-a2) - S * carry * normalCDF(-a1) +
               S * discR * coef *
                   (-normalCDF(a1 - reflT) + ebT * normalCDF(a1));
    }

    // Fixed strike: strike against the realised extreme.
    const double d1 = (std::log(S / K) + (b + 0.5 * sigma * sigma) * T) / sigT;
    const double d2 = d1 - sigT;
    const double reflPow = std::pow(S / K, -2.0 * b / (sigma * sigma));

    if (o.isCall()) {
        // Pays max(max(S) - K, 0); the running max starts at S.
        if (K > S) {
            return S * carry * normalCDF(d1) - K * discR * normalCDF(d2) +
                   S * discR * coef *
                       (-reflPow * normalCDF(d1 - reflT) + ebT * normalCDF(d1));
        }
        // K <= S: the payoff already contains (S - K) for certain.
        return discR * (S - K) + S * carry * normalCDF(a1) - S * discR * normalCDF(a2) +
               S * discR * coef *
                   (-normalCDF(a1 - reflT) + ebT * normalCDF(a1));
    }

    // Pays max(K - min(S), 0); the running min starts at S.
    if (K < S) {
        return -S * carry * normalCDF(-d1) + K * discR * normalCDF(-d2) +
               S * discR * coef *
                   (reflPow * normalCDF(-d1 + reflT) - ebT * normalCDF(-d1));
    }
    // K >= S: the payoff already contains (K - S) for certain.
    return discR * (K - S) - S * carry * normalCDF(-a1) + S * discR * normalCDF(-a2) +
           S * discR * coef *
               (normalCDF(-a1 + reflT) - ebT * normalCDF(-a1));
}

LookbackMonteCarlo::LookbackMonteCarlo(size_t numPaths, size_t numSteps, unsigned seed)
    : numPaths_(numPaths), numSteps_(numSteps), seed_(seed) {
    if (numPaths_ == 0) throw std::invalid_argument("Lookback MC needs at least one path");
    if (numSteps_ == 0) throw std::invalid_argument("Lookback MC needs at least one step");
}

LookbackMonteCarlo::Result LookbackMonteCarlo::price(const LookbackParams& params) const {
    const OptionParams& o = params.option;
    const double dt = o.timeToMaturity / static_cast<double>(numSteps_);
    const double drift = (o.riskFreeRate - o.dividendYield - 0.5 * o.volatility * o.volatility) * dt;
    const double vol = o.volatility * std::sqrt(dt);
    const double sigma2dt = o.volatility * o.volatility * dt;
    const double discount = std::exp(-o.riskFreeRate * o.timeToMaturity);

    // Which extreme the contract watches: minima for floating calls and
    // fixed puts, maxima for floating puts and fixed calls.
    const bool tracksMin =
        (params.lookbackType == LookbackType::FLOATING_STRIKE) == o.isCall();

    std::mt19937_64 rng(seed_);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double sum = 0.0, sumSq = 0.0;
    for (size_t path = 0; path < numPaths_; ++path) {
        double logS = std::log(o.spotPrice);
        double extreme = logS;
        for (size_t step = 0; step < numSteps_; ++step) {
            const double logNext = logS + drift + vol * normal(rng);
            // The extreme of the Brownian bridge between two known
            // endpoints has an exact sampleable distribution:
            //   min = (x0 + x1 - sqrt((x1-x0)^2 - 2 sigma^2 dt ln U)) / 2
            // and symmetrically for the max. Sampling it removes the
            // discrete-monitoring bias a node-only extreme would carry.
            const double u = std::fmax(uniform(rng), 1e-300);
            const double gap = std::sqrt((logNext - logS) * (logNext - logS) -
                                         2.0 * sigma2dt * std::log(u));
            if (tracksMin) {
                extreme = std::fmin(extreme, 0.5 * (logS + logNext - gap));
            } else {
                extreme = std::fmax(extreme, 0.5 * (logS + logNext + gap));
            }
            logS = logNext;
        }
        const double sT = std::exp(logS);
        const double ext = std::exp(extreme);

        double payoff;
        if (params.lookbackType == LookbackType::FLOATING_STRIKE) {
            payoff = o.isCall() ? sT - ext : ext - sT;
        } else {
            payoff = o.isCall() ? std::fmax(ext - o.strikePrice, 0.0)
                                : std::fmax(o.strikePrice - ext, 0.0);
        }
        const double discounted = discount * payoff;
        sum += discounted;
        sumSq += discounted * discounted;
    }

    const double n = static_cast<double>(numPaths_);
    const double mean = sum / n;
    const double variance = std::fmax(sumSq / n - mean * mean, 0.0);
    return {mean, std::sqrt(variance / n)};
}

} // namespace Options
