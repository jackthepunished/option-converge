#include "options/BarrierOption.h"

#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace Options {

namespace {

double normalCDF(double x) noexcept {
    constexpr double kInvSqrt2 = 0.70710678118654752440;
    return 0.5 * std::erfc(-x * kInvSqrt2);
}

double vanillaBS(const OptionParams& p) noexcept {
    const double sqrtT = std::sqrt(p.timeToMaturity);
    const double d1 = (std::log(p.spotPrice / p.strikePrice) +
                       (p.riskFreeRate - p.dividendYield + 0.5 * p.volatility * p.volatility) *
                           p.timeToMaturity) /
                      (p.volatility * sqrtT);
    const double d2 = d1 - p.volatility * sqrtT;
    const double discS = p.spotPrice * std::exp(-p.dividendYield * p.timeToMaturity);
    const double discK = p.strikePrice * std::exp(-p.riskFreeRate * p.timeToMaturity);
    return p.isCall() ? discS * normalCDF(d1) - discK * normalCDF(d2)
                      : discK * normalCDF(-d2) - discS * normalCDF(-d1);
}

double vanillaPayoff(double spot, const OptionParams& p) noexcept {
    return p.isCall() ? std::fmax(spot - p.strikePrice, 0.0)
                      : std::fmax(p.strikePrice - spot, 0.0);
}

} // namespace

// Reiner-Rubinstein knock-in prices via the standard A/B/C/D building blocks
// (Haug's parameterisation, zero rebate). phi = +1 call / -1 put selects the
// payoff; eta = +1 down / -1 up selects the barrier side. Knock-outs follow
// from in-out parity, exact at zero rebate.
double AnalyticBarrier::price(const BarrierParams& params) {
    const OptionParams& o = params.option;

    if (params.alreadyBreached()) {
        return params.isKnockIn() ? vanillaBS(o) : 0.0;
    }

    const double S = o.spotPrice;
    const double K = o.strikePrice;
    const double B = params.barrier;
    const double r = o.riskFreeRate;
    const double q = o.dividendYield;
    const double sigma = o.volatility;
    const double T = o.timeToMaturity;

    const double sqrtT = std::sqrt(T);
    const double sigSqrtT = sigma * sqrtT;
    const double mu = (r - q) / (sigma * sigma) - 0.5;
    const double discS = S * std::exp(-q * T);
    const double discK = K * std::exp(-r * T);
    const double phi = o.isCall() ? 1.0 : -1.0;
    const double eta = params.isUp() ? -1.0 : 1.0;

    const double x1 = std::log(S / K) / sigSqrtT + (1.0 + mu) * sigSqrtT;
    const double x2 = std::log(S / B) / sigSqrtT + (1.0 + mu) * sigSqrtT;
    const double y1 = std::log(B * B / (S * K)) / sigSqrtT + (1.0 + mu) * sigSqrtT;
    const double y2 = std::log(B / S) / sigSqrtT + (1.0 + mu) * sigSqrtT;
    const double powA = std::pow(B / S, 2.0 * mu);        // (B/S)^(2mu)
    const double powB = std::pow(B / S, 2.0 * mu + 2.0);  // (B/S)^(2mu+2)

    const double A = phi * (discS * normalCDF(phi * x1) - discK * normalCDF(phi * (x1 - sigSqrtT)));
    const double Bv = phi * (discS * normalCDF(phi * x2) - discK * normalCDF(phi * (x2 - sigSqrtT)));
    const double C = phi * (discS * powB * normalCDF(eta * y1) -
                            discK * powA * normalCDF(eta * (y1 - sigSqrtT)));
    const double D = phi * (discS * powB * normalCDF(eta * y2) -
                            discK * powA * normalCDF(eta * (y2 - sigSqrtT)));

    double knockIn;
    if (o.isCall()) {
        knockIn = params.isUp() ? (K >= B ? A : Bv - C + D)   // up-and-in call
                                : (K >= B ? C : A - Bv + D);  // down-and-in call
    } else {
        knockIn = params.isUp() ? (K >= B ? A - Bv + D : C)   // up-and-in put
                                : (K >= B ? Bv - C + D : A);  // down-and-in put
    }

    return params.isKnockIn() ? knockIn : vanillaBS(o) - knockIn;
}

BarrierMonteCarlo::BarrierMonteCarlo(size_t numPaths, size_t numSteps, unsigned seed)
    : numPaths_(numPaths), numSteps_(numSteps), seed_(seed) {
    if (numPaths_ == 0) throw std::invalid_argument("Barrier MC needs at least one path");
    if (numSteps_ == 0) throw std::invalid_argument("Barrier MC needs at least one step");
}

BarrierMonteCarlo::Result BarrierMonteCarlo::price(const BarrierParams& params) const {
    const OptionParams& o = params.option;

    if (params.alreadyBreached()) {
        // Degenerate cases carry no sampling error: a breached knock-out is
        // worthless, a breached knock-in is the vanilla option (priced by MC
        // for consistency of the estimator with the non-degenerate cases).
        if (!params.isKnockIn()) return {0.0, 0.0};
    }

    const double S0 = o.spotPrice;
    const double B = params.barrier;
    const double r = o.riskFreeRate;
    const double q = o.dividendYield;
    const double sigma = o.volatility;
    const double T = o.timeToMaturity;
    const double dt = T / static_cast<double>(numSteps_);
    const double drift = (r - q - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);
    const double discount = std::exp(-r * T);
    const double sigma2dt = sigma * sigma * dt;
    const bool up = params.isUp();
    const bool breachedAtStart = params.alreadyBreached();

    std::mt19937_64 rng(seed_);
    std::normal_distribution<double> normal(0.0, 1.0);

    double sum = 0.0;
    double sumSq = 0.0;
    for (size_t path = 0; path < numPaths_; ++path) {
        double s = S0;
        // Probability the path has NOT touched the barrier, conditioned on
        // the simulated nodes. Between nodes a Brownian bridge crosses a
        // barrier both endpoints avoid with probability
        // exp(-2 ln(S_i/B) ln(S_{i+1}/B) / (sigma^2 dt)), which is what makes
        // coarse time grids price the continuously monitored contract.
        double survival = breachedAtStart ? 0.0 : 1.0;
        for (size_t step = 0; step < numSteps_; ++step) {
            const double sNext = s * std::exp(drift + vol * normal(rng));
            if (survival > 0.0) {
                const bool touched = up ? sNext >= B : sNext <= B;
                if (touched) {
                    survival = 0.0;
                } else {
                    const double crossProb =
                        std::exp(-2.0 * std::log(s / B) * std::log(sNext / B) / sigma2dt);
                    survival *= 1.0 - crossProb;
                }
            }
            s = sNext;
        }
        const double weight = params.isKnockIn() ? 1.0 - survival : survival;
        const double payoff = discount * weight * vanillaPayoff(s, o);
        sum += payoff;
        sumSq += payoff * payoff;
    }

    const double n = static_cast<double>(numPaths_);
    const double mean = sum / n;
    const double variance = std::fmax(sumSq / n - mean * mean, 0.0);
    return {mean, std::sqrt(variance / n)};
}

} // namespace Options
