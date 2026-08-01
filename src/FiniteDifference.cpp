#include "options/FiniteDifference.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Options {

FiniteDifference::FiniteDifference(size_t spotSteps, size_t timeSteps, FDScheme scheme)
    : spotSteps_(spotSteps), timeSteps_(timeSteps), scheme_(scheme) {
    if (spotSteps_ < 4) {
        throw std::invalid_argument("Finite difference grid needs at least 4 spot steps");
    }
    if (timeSteps_ < 1) {
        throw std::invalid_argument("Finite difference grid needs at least 1 time step");
    }
    // An even interval count puts the initial spot exactly on a grid node, so
    // the price is read off directly instead of interpolated.
    if (spotSteps_ % 2 != 0) {
        ++spotSteps_;
    }
}

std::string FiniteDifference::getName() const {
    const char* scheme = scheme_ == FDScheme::EXPLICIT   ? "explicit"
                         : scheme_ == FDScheme::IMPLICIT ? "implicit"
                                                         : "Crank-Nicolson";
    return "Finite Difference (" + std::string(scheme) + ", " +
           std::to_string(spotSteps_) + "x" + std::to_string(timeSteps_) + ")";
}

PricingResult FiniteDifference::price(const OptionParams& params) {
    if (params.hasDiscreteDividends()) {
        // European contracts take the escrowed-dividend transformation; the
        // American grid would need a time-dependent exercise payoff, which
        // the lattice already provides - refuse rather than approximate.
        if (params.isAmerican()) {
            throw std::invalid_argument(
                "Finite difference does not price American options with discrete "
                "dividends; use BinomialTree.");
        }
        return price(params.escrowed());
    }

    const auto start = std::chrono::high_resolution_clock::now();

    PricingResult result;
    result.price = solveGrid(params, result.memoryUsed);

    const auto end = std::chrono::high_resolution_clock::now();
    result.computationTime =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    return result;
}

double FiniteDifference::payoff(double spot, const OptionParams& params) noexcept {
    return params.isCall() ? std::fmax(spot - params.strikePrice, 0.0)
                           : std::fmax(params.strikePrice - spot, 0.0);
}

double FiniteDifference::lowerBoundary(double spot, double tau,
                                       const OptionParams& params) noexcept {
    if (params.isCall()) {
        return 0.0;
    }
    if (params.isAmerican()) {
        // Deep in the money the American put is exercised immediately.
        return params.strikePrice - spot;
    }
    return params.strikePrice * std::exp(-params.riskFreeRate * tau) -
           spot * std::exp(-params.dividendYield * tau);
}

double FiniteDifference::upperBoundary(double spot, double tau,
                                       const OptionParams& params) noexcept {
    if (params.isPut()) {
        return 0.0;
    }
    if (params.isAmerican()) {
        return std::fmax(spot - params.strikePrice, 0.0);
    }
    return spot * std::exp(-params.dividendYield * tau) -
           params.strikePrice * std::exp(-params.riskFreeRate * tau);
}

double FiniteDifference::solveGrid(const OptionParams& params, size_t& memoryUsed) const {
    const double S = params.spotPrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;

    // Log-spot grid centred on ln(S) so the spot sits on the middle node.
    // Width covers 5 standard deviations past the strike's log-distance,
    // keeping both the payoff kink and the tails inside the domain.
    const size_t M = spotSteps_;
    const double halfWidth =
        std::fabs(std::log(params.strikePrice / S)) + 5.0 * sigma * std::sqrt(T) + 1e-8;
    const double x0 = std::log(S);
    const double dx = 2.0 * halfWidth / static_cast<double>(M);

    // The explicit scheme is only stable while the discrete update keeps
    // positive coefficients: dt <= 1 / (sigma^2/dx^2 + r). Refine the time
    // axis to honour that rather than returning an oscillating solution.
    size_t N = timeSteps_;
    if (scheme_ == FDScheme::EXPLICIT) {
        const double dtMax = 1.0 / (sigma * sigma / (dx * dx) + r);
        const auto needed = static_cast<size_t>(std::ceil(T / dtMax));
        N = std::max(N, needed + 1);
    }
    const double dt = T / static_cast<double>(N);

    // Constant theta-scheme stencil: dt * L where
    // L V = a V_xx + b V_x - r V,  a = sigma^2/2,  b = r - q - sigma^2/2.
    const double a = 0.5 * sigma * sigma;
    const double b = r - q - a;
    const double lo = dt * (a / (dx * dx) - b / (2.0 * dx));
    const double di = dt * (-2.0 * a / (dx * dx) - r);
    const double up = dt * (a / (dx * dx) + b / (2.0 * dx));
    const double theta = scheme_ == FDScheme::EXPLICIT   ? 0.0
                         : scheme_ == FDScheme::IMPLICIT ? 1.0
                                                         : 0.5;

    std::vector<double> spot(M + 1);
    for (size_t j = 0; j <= M; ++j) {
        spot[j] = std::exp(x0 - halfWidth + static_cast<double>(j) * dx);
    }

    std::vector<double> values(M + 1);
    for (size_t j = 0; j <= M; ++j) {
        values[j] = payoff(spot[j], params);
    }

    std::vector<double> rhs(M + 1);
    std::vector<double> scratch(M + 1);  // Thomas forward-sweep workspace
    memoryUsed = 4 * (M + 1) * sizeof(double);

    // March backward in time: level n has time-to-maturity tau = T - n*dt.
    for (size_t n = N; n-- > 0;) {
        const double tau = T - static_cast<double>(n) * dt;

        for (size_t j = 1; j < M; ++j) {
            rhs[j] = values[j] + (1.0 - theta) * (lo * values[j - 1] + di * values[j] +
                                                  up * values[j + 1]);
        }
        const double lowerVal = lowerBoundary(spot[0], tau, params);
        const double upperVal = upperBoundary(spot[M], tau, params);

        if (theta == 0.0) {
            for (size_t j = 1; j < M; ++j) {
                values[j] = rhs[j];
            }
        } else {
            // Thomas algorithm on the interior nodes. The matrix
            // (I - theta*dt*L) is constant and diagonally dominant, so the
            // sweep is unconditionally well-posed. Dirichlet boundaries fold
            // into the first and last right-hand sides.
            const double subDiag = -theta * lo;
            const double diag = 1.0 - theta * di;
            const double superDiag = -theta * up;
            rhs[1] -= subDiag * lowerVal;
            rhs[M - 1] -= superDiag * upperVal;

            scratch[1] = superDiag / diag;
            rhs[1] /= diag;
            for (size_t j = 2; j < M; ++j) {
                const double denom = diag - subDiag * scratch[j - 1];
                scratch[j] = superDiag / denom;
                rhs[j] = (rhs[j] - subDiag * rhs[j - 1]) / denom;
            }
            values[M - 1] = rhs[M - 1];
            for (size_t j = M - 1; j-- > 1;) {
                values[j] = rhs[j] - scratch[j] * values[j + 1];
            }
        }

        values[0] = lowerVal;
        values[M] = upperVal;

        // Early-exercise projection: at every node the holder takes the
        // better of continuing and exercising now.
        if (params.isAmerican()) {
            for (size_t j = 0; j <= M; ++j) {
                values[j] = std::fmax(values[j], payoff(spot[j], params));
            }
        }
    }

    return values[M / 2];  // The spot sits exactly on the centre node.
}

} // namespace Options
