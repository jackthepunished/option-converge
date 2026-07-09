#include "options/MonteCarlo.h"
#include "options/BlackScholes.h"

#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace Options {

namespace {

// Time steps per simulated path. Vanilla European payoffs only depend on the
// terminal price, but the discretization schemes are exercised as specified so
// their bias can be studied by the convergence analyzer.
constexpr size_t kStepsPerPath = 252;

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

[[nodiscard]] double mean(const std::vector<double>& xs) {
    return std::accumulate(xs.begin(), xs.end(), 0.0) / static_cast<double>(xs.size());
}

} // namespace

MonteCarlo::MonteCarlo(size_t numPaths, DiscretizationScheme scheme,
                       VarianceReduction varRed, unsigned int seed)
    : numPaths_(numPaths), scheme_(scheme), varRed_(varRed), rng_(seed) {
    if (numPaths == 0) {
        throw std::invalid_argument("Monte Carlo requires at least one path");
    }
}

void MonteCarlo::setSeed(unsigned int seed) {
    rng_.seed(seed);
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

double MonteCarlo::randomNormal() {
    std::normal_distribution<double> normal(0.0, 1.0);
    return normal(rng_);
}

double MonteCarlo::payoff(double finalPrice, const OptionParams& params) const noexcept {
    return params.isCall() ? std::fmax(finalPrice - params.strikePrice, 0.0)
                           : std::fmax(params.strikePrice - finalPrice, 0.0);
}

double MonteCarlo::simulatePathEuler(const OptionParams& params) {
    std::vector<double> normals(kStepsPerPath);
    for (double& z : normals) {
        z = randomNormal();
    }
    return terminalPrice(params, normals, DiscretizationScheme::EULER, false);
}

double MonteCarlo::simulatePathMilstein(const OptionParams& params) {
    std::vector<double> normals(kStepsPerPath);
    for (double& z : normals) {
        z = randomNormal();
    }
    return terminalPrice(params, normals, DiscretizationScheme::MILSTEIN, false);
}

double MonteCarlo::calculateStandardError(const std::vector<double>& payoffs) const {
    const size_t n = payoffs.size();
    if (n < 2) {
        return 0.0;
    }
    const double m = mean(payoffs);
    double sumSq = 0.0;
    for (double x : payoffs) {
        const double dev = x - m;
        sumSq += dev * dev;
    }
    const double variance = sumSq / static_cast<double>(n - 1);
    return std::sqrt(variance / static_cast<double>(n));
}

// Plain Monte Carlo: one independent path per sample.
double MonteCarlo::priceBasic(const OptionParams& params, double& stdError) {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);

    std::vector<double> discounted(numPaths_);
    for (size_t i = 0; i < numPaths_; ++i) {
        const double sT = (scheme_ == DiscretizationScheme::EULER)
                              ? simulatePathEuler(params)
                              : simulatePathMilstein(params);
        discounted[i] = discount * payoff(sT, params);
    }

    stdError = calculateStandardError(discounted);
    return mean(discounted);
}

// Antithetic variates: each draw is reused with its sign flipped, and the two
// path payoffs are averaged. The pair averages are the i.i.d. samples.
double MonteCarlo::priceAntithetic(const OptionParams& params, double& stdError) {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);
    const size_t numPairs = (numPaths_ + 1) / 2;

    std::vector<double> normals(kStepsPerPath);
    std::vector<double> pairMeans(numPairs);

    for (size_t i = 0; i < numPairs; ++i) {
        for (double& z : normals) {
            z = randomNormal();
        }
        const double sUp = terminalPrice(params, normals, scheme_, false);
        const double sDown = terminalPrice(params, normals, scheme_, true);
        pairMeans[i] = discount * 0.5 * (payoff(sUp, params) + payoff(sDown, params));
    }

    stdError = calculateStandardError(pairMeans);
    return mean(pairMeans);
}

// Control variates: the discounted terminal stock price has known expectation
// S0 * exp(-q T), so its sampling error is subtracted from the payoff estimate
// with the variance-minimizing coefficient beta = Cov(Y, C) / Var(C).
double MonteCarlo::priceControlVariate(const OptionParams& params, double& stdError) {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);
    const double controlMean =
        params.spotPrice * std::exp(-params.dividendYield * params.timeToMaturity);

    std::vector<double> discounted(numPaths_);
    std::vector<double> controls(numPaths_);
    for (size_t i = 0; i < numPaths_; ++i) {
        const double sT = (scheme_ == DiscretizationScheme::EULER)
                              ? simulatePathEuler(params)
                              : simulatePathMilstein(params);
        discounted[i] = discount * payoff(sT, params);
        controls[i] = discount * sT;
    }

    const double yBar = mean(discounted);
    const double cBar = mean(controls);
    double covYC = 0.0;
    double varC = 0.0;
    for (size_t i = 0; i < numPaths_; ++i) {
        covYC += (discounted[i] - yBar) * (controls[i] - cBar);
        varC += (controls[i] - cBar) * (controls[i] - cBar);
    }
    const double beta = (varC > 0.0) ? covYC / varC : 0.0;

    std::vector<double> adjusted(numPaths_);
    for (size_t i = 0; i < numPaths_; ++i) {
        adjusted[i] = discounted[i] - beta * (controls[i] - controlMean);
    }

    stdError = calculateStandardError(adjusted);
    return mean(adjusted);
}

// Antithetic pairs first, then the control variate applied to the pair means.
double MonteCarlo::priceBoth(const OptionParams& params, double& stdError) {
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);
    const double controlMean =
        params.spotPrice * std::exp(-params.dividendYield * params.timeToMaturity);
    const size_t numPairs = (numPaths_ + 1) / 2;

    std::vector<double> normals(kStepsPerPath);
    std::vector<double> pairPayoffs(numPairs);
    std::vector<double> pairControls(numPairs);

    for (size_t i = 0; i < numPairs; ++i) {
        for (double& z : normals) {
            z = randomNormal();
        }
        const double sUp = terminalPrice(params, normals, scheme_, false);
        const double sDown = terminalPrice(params, normals, scheme_, true);
        pairPayoffs[i] = discount * 0.5 * (payoff(sUp, params) + payoff(sDown, params));
        pairControls[i] = discount * 0.5 * (sUp + sDown);
    }

    const double yBar = mean(pairPayoffs);
    const double cBar = mean(pairControls);
    double covYC = 0.0;
    double varC = 0.0;
    for (size_t i = 0; i < numPairs; ++i) {
        covYC += (pairPayoffs[i] - yBar) * (pairControls[i] - cBar);
        varC += (pairControls[i] - cBar) * (pairControls[i] - cBar);
    }
    const double beta = (varC > 0.0) ? covYC / varC : 0.0;

    std::vector<double> adjusted(numPairs);
    for (size_t i = 0; i < numPairs; ++i) {
        adjusted[i] = pairPayoffs[i] - beta * (pairControls[i] - controlMean);
    }

    stdError = calculateStandardError(adjusted);
    return mean(adjusted);
}

PricingResult MonteCarlo::price(const OptionParams& params) {
    if (!params.isEuropean()) {
        throw std::invalid_argument(
            "Monte Carlo engine only prices European options. "
            "Use Binomial Tree for American options.");
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

    // Sample vectors plus the per-path normal draws dominate.
    result.memoryUsed = numPaths_ * 2 * sizeof(double) + kStepsPerPath * sizeof(double);

    return result;
}

Greeks MonteCarlo::calculateGreeks(const OptionParams& params) {
    // Finite differences on independent samples would be dominated by Monte
    // Carlo noise. Restoring the same RNG state before each repricing (common
    // random numbers) makes the noise cancel in the differences.
    const std::mt19937 saved = rng_;
    const auto priceAt = [&](const OptionParams& p) {
        rng_ = saved;
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

    rng_ = saved;
    return greeks;
}

} // namespace Options
