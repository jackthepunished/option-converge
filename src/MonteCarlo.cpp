#include "options/MonteCarlo.h"

#include <chrono>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace Options {

namespace {

// Time steps per simulated path. Vanilla European payoffs only depend on the
// terminal price, but the discretization schemes are exercised as specified so
// their bias can be studied by the convergence analyzer.
constexpr size_t kStepsPerPath = 252;

// Path work is split into this many chunks regardless of how many threads
// run them. Each chunk owns a deterministically seeded RNG, so the estimate
// depends only on the seed and the partition - never on the thread count.
constexpr int kNumChunks = 256;

// Distinct, deterministic seed per (engine seed, chunk) pair. The golden
// ratio increment keeps neighbouring chunk seeds far apart in state space.
[[nodiscard]] unsigned int chunkSeed(unsigned int seed, int chunk) noexcept {
    return seed + 0x9E3779B9u * static_cast<unsigned int>(chunk + 1);
}

// Number of samples chunk c takes when n samples are spread over the chunks.
[[nodiscard]] size_t chunkShare(size_t n, int c) noexcept {
    const size_t base = n / kNumChunks;
    const size_t rem = n % kNumChunks;
    return base + (static_cast<size_t>(c) < rem ? 1 : 0);
}

// Advance one GBM path to maturity given a pre-drawn vector of standard
// normals. Milstein adds the second-order correction term to Euler.
[[nodiscard]] double terminalPrice(const OptionParams& params,
                                   const std::vector<double>& normals,
                                   DiscretizationScheme scheme,
                                   bool antithetic) {
    const double dt = params.timeToMaturity / static_cast<double>(normals.size());
    const double sqrtDt = std::sqrt(dt);
    const double drift = params.riskFreeRate - params.dividendYield;
    const double sigma = params.volatility;

    double s = params.spotPrice;
    for (double z : normals) {
        if (antithetic) {
            z = -z;
        }
        const double dW = sqrtDt * z;
        double ds = s * (drift * dt + sigma * dW);
        if (scheme == DiscretizationScheme::MILSTEIN) {
            ds += 0.5 * sigma * sigma * s * (dW * dW - dt);
        }
        s += ds;
        // Euler can step through zero for large negative draws; GBM never
        // goes negative, so clamp to keep the path meaningful.
        s = std::fmax(s, 1e-12);
    }
    return s;
}

// Sample standard error from accumulated first and second moments (n-1
// denominator, matching the previous vector-based implementation).
[[nodiscard]] double standardErrorFromSums(double sum, double sumSq, size_t n) noexcept {
    if (n < 2) {
        return 0.0;
    }
    const double nd = static_cast<double>(n);
    const double mean = sum / nd;
    const double variance = std::fmax(sumSq - nd * mean * mean, 0.0) / (nd - 1.0);
    return std::sqrt(variance / nd);
}

} // namespace

MonteCarlo::MonteCarlo(size_t numPaths, DiscretizationScheme scheme,
                       VarianceReduction varRed, unsigned int seed)
    : numPaths_(numPaths), scheme_(scheme), varRed_(varRed), seed_(seed) {
    if (numPaths == 0) {
        throw std::invalid_argument("Monte Carlo requires at least one path");
    }
}

void MonteCarlo::setSeed(unsigned int seed) {
    seed_ = seed;
}

std::string MonteCarlo::getName() const {
    std::string name = "Monte Carlo (" + std::to_string(numPaths_) + " paths, ";
    name += (scheme_ == DiscretizationScheme::EULER) ? "Euler" : "Milstein";
    switch (varRed_) {
        case VarianceReduction::NONE: name += ")"; break;
        case VarianceReduction::ANTITHETIC: name += ", antithetic)"; break;
        case VarianceReduction::CONTROL_VARIATE: name += ", control variate)"; break;
        case VarianceReduction::BOTH: name += ", antithetic + control variate)"; break;
    }
    return name;
}

double MonteCarlo::payoff(double finalPrice, const OptionParams& params) const noexcept {
    return params.isCall() ? std::fmax(finalPrice - params.strikePrice, 0.0)
                           : std::fmax(params.strikePrice - finalPrice, 0.0);
}

// Plain Monte Carlo: one independent path per sample.
double MonteCarlo::priceBasic(const OptionParams& params, double& stdError) const {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);

    double sumY = 0.0, sumY2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sumY, sumY2) schedule(static)
#endif
    for (int c = 0; c < kNumChunks; ++c) {
        std::mt19937 rng(chunkSeed(seed_, c));
        std::normal_distribution<double> normal(0.0, 1.0);
        std::vector<double> normals(kStepsPerPath);
        const size_t share = chunkShare(numPaths_, c);
        for (size_t i = 0; i < share; ++i) {
            for (double& z : normals) {
                z = normal(rng);
            }
            const double sT = terminalPrice(params, normals, scheme_, false);
            const double y = discount * payoff(sT, params);
            sumY += y;
            sumY2 += y * y;
        }
    }

    stdError = standardErrorFromSums(sumY, sumY2, numPaths_);
    return sumY / static_cast<double>(numPaths_);
}

// Antithetic variates: each draw is reused with its sign flipped, and the two
// path payoffs are averaged. The pair averages are the i.i.d. samples.
double MonteCarlo::priceAntithetic(const OptionParams& params, double& stdError) const {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);
    const size_t numPairs = (numPaths_ + 1) / 2;

    double sumY = 0.0, sumY2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sumY, sumY2) schedule(static)
#endif
    for (int c = 0; c < kNumChunks; ++c) {
        std::mt19937 rng(chunkSeed(seed_, c));
        std::normal_distribution<double> normal(0.0, 1.0);
        std::vector<double> normals(kStepsPerPath);
        const size_t share = chunkShare(numPairs, c);
        for (size_t i = 0; i < share; ++i) {
            for (double& z : normals) {
                z = normal(rng);
            }
            const double sUp = terminalPrice(params, normals, scheme_, false);
            const double sDown = terminalPrice(params, normals, scheme_, true);
            const double y = discount * 0.5 * (payoff(sUp, params) + payoff(sDown, params));
            sumY += y;
            sumY2 += y * y;
        }
    }

    stdError = standardErrorFromSums(sumY, sumY2, numPairs);
    return sumY / static_cast<double>(numPairs);
}

// Control variates: the discounted terminal stock price has known expectation
// S0 * exp(-q T), so its sampling error is subtracted from the payoff estimate
// with the variance-minimizing coefficient beta = Cov(Y, C) / Var(C). All five
// moment sums are additive, so the correction is applied after the reduction.
double MonteCarlo::priceControlVariate(const OptionParams& params, double& stdError) const {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);
    const double controlMean =
        params.spotPrice * std::exp(-params.dividendYield * params.timeToMaturity);

    double sumY = 0.0, sumY2 = 0.0, sumC = 0.0, sumC2 = 0.0, sumYC = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sumY, sumY2, sumC, sumC2, sumYC) schedule(static)
#endif
    for (int c = 0; c < kNumChunks; ++c) {
        std::mt19937 rng(chunkSeed(seed_, c));
        std::normal_distribution<double> normal(0.0, 1.0);
        std::vector<double> normals(kStepsPerPath);
        const size_t share = chunkShare(numPaths_, c);
        for (size_t i = 0; i < share; ++i) {
            for (double& z : normals) {
                z = normal(rng);
            }
            const double sT = terminalPrice(params, normals, scheme_, false);
            const double y = discount * payoff(sT, params);
            const double ctrl = discount * sT;
            sumY += y;
            sumY2 += y * y;
            sumC += ctrl;
            sumC2 += ctrl * ctrl;
            sumYC += y * ctrl;
        }
    }

    const double n = static_cast<double>(numPaths_);
    const double meanY = sumY / n;
    const double meanC = sumC / n;
    const double covYC = sumYC - n * meanY * meanC;
    const double varC = std::fmax(sumC2 - n * meanC * meanC, 0.0);
    const double beta = (varC > 0.0) ? covYC / varC : 0.0;

    const double varY = std::fmax(sumY2 - n * meanY * meanY, 0.0);
    const double varAdj =
        std::fmax(varY - 2.0 * beta * covYC + beta * beta * varC, 0.0) / (n - 1.0);
    stdError = std::sqrt(varAdj / n);
    return meanY - beta * (meanC - controlMean);
}

// Antithetic pairs first, then the control variate applied to the pair means.
double MonteCarlo::priceBoth(const OptionParams& params, double& stdError) const {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);
    const double controlMean =
        params.spotPrice * std::exp(-params.dividendYield * params.timeToMaturity);
    const size_t numPairs = (numPaths_ + 1) / 2;

    double sumY = 0.0, sumY2 = 0.0, sumC = 0.0, sumC2 = 0.0, sumYC = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sumY, sumY2, sumC, sumC2, sumYC) schedule(static)
#endif
    for (int c = 0; c < kNumChunks; ++c) {
        std::mt19937 rng(chunkSeed(seed_, c));
        std::normal_distribution<double> normal(0.0, 1.0);
        std::vector<double> normals(kStepsPerPath);
        const size_t share = chunkShare(numPairs, c);
        for (size_t i = 0; i < share; ++i) {
            for (double& z : normals) {
                z = normal(rng);
            }
            const double sUp = terminalPrice(params, normals, scheme_, false);
            const double sDown = terminalPrice(params, normals, scheme_, true);
            const double y = discount * 0.5 * (payoff(sUp, params) + payoff(sDown, params));
            const double ctrl = discount * 0.5 * (sUp + sDown);
            sumY += y;
            sumY2 += y * y;
            sumC += ctrl;
            sumC2 += ctrl * ctrl;
            sumYC += y * ctrl;
        }
    }

    const double n = static_cast<double>(numPairs);
    const double meanY = sumY / n;
    const double meanC = sumC / n;
    const double covYC = sumYC - n * meanY * meanC;
    const double varC = std::fmax(sumC2 - n * meanC * meanC, 0.0);
    const double beta = (varC > 0.0) ? covYC / varC : 0.0;

    const double varY = std::fmax(sumY2 - n * meanY * meanY, 0.0);
    const double varAdj =
        std::fmax(varY - 2.0 * beta * covYC + beta * beta * varC, 0.0) / (n - 1.0);
    stdError = std::sqrt(varAdj / n);
    return meanY - beta * (meanC - controlMean);
}

PricingResult MonteCarlo::price(const OptionParams& params) {
    if (!params.isEuropean()) {
        throw std::invalid_argument(
            "Monte Carlo engine only prices European options. "
            "Use Binomial Tree for American options.");
    }

    // Escrowed-dividend model: simulate the dividend-stripped spot.
    if (params.hasDiscreteDividends()) {
        return price(params.escrowed());
    }

    const auto start = std::chrono::high_resolution_clock::now();

    PricingResult result;
    switch (varRed_) {
        case VarianceReduction::NONE:
            result.price = priceBasic(params, result.standardError);
            break;
        case VarianceReduction::ANTITHETIC:
            result.price = priceAntithetic(params, result.standardError);
            break;
        case VarianceReduction::CONTROL_VARIATE:
            result.price = priceControlVariate(params, result.standardError);
            break;
        case VarianceReduction::BOTH:
            result.price = priceBoth(params, result.standardError);
            break;
    }

    const auto end = std::chrono::high_resolution_clock::now();
    result.computationTime =
        std::chrono::duration<double, std::milli>(end - start).count();

    // Accumulation happens in scalar sums; the per-chunk normal draw buffers
    // dominate what little the engine allocates.
    result.memoryUsed = kNumChunks * kStepsPerPath * sizeof(double);

    return result;
}

Greeks MonteCarlo::calculateGreeks(const OptionParams& params) {
    // Common random numbers by construction: pricing is a pure function of
    // the seed, so every repricing below reuses exactly the same draws and
    // the sampling noise cancels in the finite differences.
    const auto priceAt = [&](const OptionParams& p) {
        return price(p).price;
    };

    const double hS = 0.01;
    const double hSigma = 0.01;
    const double hR = 0.01;
    const double hT = 1.0 / 365.0;

    Greeks greeks;

    OptionParams up = params;
    OptionParams down = params;

    up.spotPrice = params.spotPrice + hS;
    down.spotPrice = params.spotPrice - hS;
    const double priceUp = priceAt(up);
    const double priceDown = priceAt(down);
    const double priceCenter = priceAt(params);
    greeks.delta = (priceUp - priceDown) / (2.0 * hS);
    greeks.gamma = (priceUp - 2.0 * priceCenter + priceDown) / (hS * hS);

    // Vega/rho per 1% and theta per day, matching the analytic engine.
    up = params;
    down = params;
    up.volatility = params.volatility + hSigma;
    down.volatility = params.volatility - hSigma;
    greeks.vega = (priceAt(up) - priceAt(down)) / (2.0 * hSigma) / 100.0;

    up = params;
    down = params;
    up.riskFreeRate = params.riskFreeRate + hR;
    down.riskFreeRate = params.riskFreeRate - hR;
    greeks.rho = (priceAt(up) - priceAt(down)) / (2.0 * hR) / 100.0;

    if (params.timeToMaturity > hT) {
        OptionParams forward = params;
        forward.timeToMaturity = params.timeToMaturity - hT;
        greeks.theta = (priceAt(forward) - priceCenter) / hT / 365.0;
    }

    return greeks;
}

} // namespace Options
