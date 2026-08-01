#include "options/AsianOption.h"

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

double vanillaPayoff(double average, const OptionParams& p) noexcept {
    return p.isCall() ? std::fmax(average - p.strikePrice, 0.0)
                      : std::fmax(p.strikePrice - average, 0.0);
}

// Log-moments of the discrete geometric average with n uniform fixings.
// Sum_{i,j} min(i,j) = n(n+1)(2n+1)/6 gives the variance term.
void geometricLogMoments(const OptionParams& o, size_t n, double& mean, double& variance) noexcept {
    const double nd = static_cast<double>(n);
    const double T = o.timeToMaturity;
    mean = std::log(o.spotPrice) +
           (o.riskFreeRate - o.dividendYield - 0.5 * o.volatility * o.volatility) * T *
               (nd + 1.0) / (2.0 * nd);
    variance = o.volatility * o.volatility * T * (nd + 1.0) * (2.0 * nd + 1.0) / (6.0 * nd * nd);
}

} // namespace

double AnalyticGeometricAsian::price(const AsianParams& params) {
    if (params.averaging != AveragingType::GEOMETRIC) {
        throw std::invalid_argument(
            "No closed form exists for the arithmetic average; use AsianMonteCarlo.");
    }
    const OptionParams& o = params.option;

    double m, v2;
    geometricLogMoments(o, params.numFixings, m, v2);
    const double v = std::sqrt(v2);
    const double discount = std::exp(-o.riskFreeRate * o.timeToMaturity);
    const double forwardG = std::exp(m + 0.5 * v2);  // E[G] under Q
    const double d1 = (m - std::log(o.strikePrice) + v2) / v;
    const double d2 = d1 - v;

    return o.isCall()
               ? discount * (forwardG * normalCDF(d1) - o.strikePrice * normalCDF(d2))
               : discount * (o.strikePrice * normalCDF(-d2) - forwardG * normalCDF(-d1));
}

AsianMonteCarlo::AsianMonteCarlo(size_t numPaths, unsigned seed)
    : numPaths_(numPaths), seed_(seed) {
    if (numPaths_ == 0) throw std::invalid_argument("Asian MC needs at least one path");
}

AsianMonteCarlo::Result AsianMonteCarlo::price(const AsianParams& params,
                                               bool useControlVariate) const {
    const OptionParams& o = params.option;
    const size_t n = params.numFixings;
    const double dt = o.timeToMaturity / static_cast<double>(n);
    const double drift = (o.riskFreeRate - o.dividendYield - 0.5 * o.volatility * o.volatility) * dt;
    const double vol = o.volatility * std::sqrt(dt);
    const double discount = std::exp(-o.riskFreeRate * o.timeToMaturity);
    const bool arithmetic = params.averaging == AveragingType::ARITHMETIC;

    std::mt19937_64 rng(seed_);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Y is the estimator target; C is the geometric control alongside it.
    double sumY = 0.0, sumY2 = 0.0;
    double sumC = 0.0, sumC2 = 0.0, sumYC = 0.0;

    for (size_t path = 0; path < numPaths_; ++path) {
        double logS = std::log(o.spotPrice);
        double sumSpot = 0.0;
        double sumLog = 0.0;
        for (size_t i = 0; i < n; ++i) {
            logS += drift + vol * normal(rng);
            sumSpot += std::exp(logS);
            sumLog += logS;
        }
        const double arithMean = sumSpot / static_cast<double>(n);
        const double geoMean = std::exp(sumLog / static_cast<double>(n));

        const double y = discount * vanillaPayoff(arithmetic ? arithMean : geoMean, o);
        sumY += y;
        sumY2 += y * y;
        if (arithmetic && useControlVariate) {
            const double c = discount * vanillaPayoff(geoMean, o);
            sumC += c;
            sumC2 += c * c;
            sumYC += y * c;
        }
    }

    const double N = static_cast<double>(numPaths_);
    const double meanY = sumY / N;
    const double varY = std::fmax(sumY2 / N - meanY * meanY, 0.0);

    if (!(arithmetic && useControlVariate)) {
        return {meanY, std::sqrt(varY / N)};
    }

    // Control-variate correction against the exact discrete geometric price.
    // Beta is estimated in-sample, the same convention the vanilla MonteCarlo
    // engine uses (small O(1/N) bias, documented in the README).
    const double meanC = sumC / N;
    const double varC = std::fmax(sumC2 / N - meanC * meanC, 0.0);
    if (varC <= 0.0) {
        return {meanY, std::sqrt(varY / N)};
    }
    const double covYC = sumYC / N - meanY * meanC;
    const double beta = covYC / varC;

    AsianParams geoParams(o, AveragingType::GEOMETRIC, n);
    const double exactGeo = AnalyticGeometricAsian::price(geoParams);

    const double adjusted = meanY - beta * (meanC - exactGeo);
    const double varAdj = std::fmax(varY - 2.0 * beta * covYC + beta * beta * varC, 0.0);
    return {adjusted, std::sqrt(varAdj / N)};
}

} // namespace Options
