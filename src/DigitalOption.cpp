#include "options/DigitalOption.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace Options {

namespace {

double normalCDF(double x) noexcept {
    constexpr double kInvSqrt2 = 0.70710678118654752440;
    return 0.5 * std::erfc(-x * kInvSqrt2);
}

double normalPDF(double x) noexcept {
    constexpr double kInvSqrt2Pi = 0.3989422804014327;
    return kInvSqrt2Pi * std::exp(-0.5 * x * x);
}

struct DValues {
    double d1;
    double d2;
};

DValues dValues(const OptionParams& o) noexcept {
    const double sqrtT = std::sqrt(o.timeToMaturity);
    const double d1 =
        (std::log(o.spotPrice / o.strikePrice) +
         (o.riskFreeRate - o.dividendYield + 0.5 * o.volatility * o.volatility) *
             o.timeToMaturity) /
        (o.volatility * sqrtT);
    return {d1, d1 - o.volatility * sqrtT};
}

} // namespace

double AnalyticDigital::price(const DigitalParams& params) {
    const OptionParams& o = params.option;
    const auto [d1, d2] = dValues(o);
    const double sign = o.isCall() ? 1.0 : -1.0;

    if (params.digitalType == DigitalType::CASH_OR_NOTHING) {
        return params.cashPayout * std::exp(-o.riskFreeRate * o.timeToMaturity) *
               normalCDF(sign * d2);
    }
    return o.spotPrice * std::exp(-o.dividendYield * o.timeToMaturity) *
           normalCDF(sign * d1);
}

double AnalyticDigital::delta(const DigitalParams& params) {
    const OptionParams& o = params.option;
    const auto [d1, d2] = dValues(o);
    const double sqrtT = std::sqrt(o.timeToMaturity);
    const double sign = o.isCall() ? 1.0 : -1.0;

    if (params.digitalType == DigitalType::CASH_OR_NOTHING) {
        // d(Q e^{-rT} N(+-d2))/dS: all the mass sits in the density at d2.
        return sign * params.cashPayout * std::exp(-o.riskFreeRate * o.timeToMaturity) *
               normalPDF(d2) / (o.spotPrice * o.volatility * sqrtT);
    }
    // d(S e^{-qT} N(+-d1))/dS: the vanilla-like N term plus the density term.
    const double carry = std::exp(-o.dividendYield * o.timeToMaturity);
    return carry * (normalCDF(sign * d1) + sign * normalPDF(d1) / (o.volatility * sqrtT));
}

DigitalLattice::DigitalLattice(size_t steps) : steps_(steps) {
    if (steps_ == 0) {
        throw std::invalid_argument("Digital lattice requires at least one step");
    }
}

double DigitalLattice::price(const DigitalParams& params) const {
    const OptionParams& o = params.option;
    const size_t n = steps_;

    const double dt = o.timeToMaturity / static_cast<double>(n);
    const double u = std::exp(o.volatility * std::sqrt(dt));
    const double d = 1.0 / u;
    const double growth = std::exp((o.riskFreeRate - o.dividendYield) * dt);
    const double p = (growth - d) / (u - d);
    if (p < 0.0 || p > 1.0) {
        throw std::invalid_argument(
            "Digital lattice is not arbitrage-free with the given step count; "
            "increase the number of steps.");
    }
    const double discount = std::exp(-o.riskFreeRate * dt);

    // Terminal spots by incremental update, as in BinomialTree.
    std::vector<double> spots(n + 1);
    const double u2 = u * u;
    spots[0] = o.spotPrice * std::pow(d, static_cast<double>(n));
    for (size_t i = 1; i <= n; ++i) {
        spots[i] = spots[i - 1] * u2;
    }

    std::vector<double> values(n + 1);
    for (size_t i = 0; i <= n; ++i) {
        const bool inTheMoney = o.isCall() ? spots[i] > o.strikePrice
                                           : spots[i] < o.strikePrice;
        if (params.digitalType == DigitalType::CASH_OR_NOTHING) {
            values[i] = inTheMoney ? params.cashPayout : 0.0;
        } else {
            values[i] = inTheMoney ? spots[i] : 0.0;
        }
    }

    const double wUp = discount * p;
    const double wDown = discount * (1.0 - p);
    for (size_t step = n; step-- > 0;) {
        for (size_t i = 0; i <= step; ++i) {
            values[i] = wUp * values[i + 1] + wDown * values[i];
        }
    }
    return values[0];
}

} // namespace Options
