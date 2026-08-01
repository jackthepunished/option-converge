#include "options/Calibration.h"

#include "options/ImpliedVolatility.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Options {

namespace {

// Parameter box for the Nelder-Mead search: wide enough to hold any
// market-plausible calibration, tight enough to keep the characteristic
// function well-behaved. Candidates outside the box score infinity, which
// the simplex treats as ordinary rejection.
constexpr double kLo[5] = {1e-6, 1e-6, 1e-3, 1e-3, -0.999};
constexpr double kHi[5] = {2.0, 2.0, 20.0, 3.0, 0.999};

bool inBounds(const double p[5]) noexcept {
    for (int i = 0; i < 5; ++i) {
        if (p[i] < kLo[i] || p[i] > kHi[i]) return false;
    }
    return true;
}

} // namespace

std::vector<ImpliedVolSurface::Point> ImpliedVolSurface::build(const MarketData& market) {
    if (market.quotes.empty()) {
        throw std::invalid_argument("Cannot build a surface from an empty quote set");
    }
    ImpliedVolatility solver;
    std::vector<Point> surface;
    surface.reserve(market.quotes.size());
    for (const MarketQuote& quote : market.quotes) {
        OptionParams params(market.spot, quote.strike, market.riskFreeRate,
                            0.2 /* ignored by the solver */, quote.maturity,
                            market.dividendYield, quote.type, ExerciseType::EUROPEAN);
        surface.push_back({quote.strike, quote.maturity,
                           solver.solve(quote.price, params).impliedVol});
    }
    return surface;
}

double HestonCalibrator::objective(const double p[5], const MarketData& market) {
    if (!inBounds(p)) {
        return std::numeric_limits<double>::infinity();
    }
    const HestonParams hp(p[0], p[1], p[2], p[3], p[4]);
    const HestonModel model(hp);

    double sse = 0.0;
    for (const MarketQuote& quote : market.quotes) {
        OptionParams params(market.spot, quote.strike, market.riskFreeRate, 0.2,
                            quote.maturity, market.dividendYield, quote.type,
                            ExerciseType::EUROPEAN);
        const double err = model.price(params) - quote.price;
        sse += quote.weight * err * err;
    }
    return sse;
}

HestonCalibrator::Result HestonCalibrator::calibrate(const MarketData& market,
                                                     const HestonParams& initialGuess) const {
    if (market.quotes.empty()) {
        throw std::invalid_argument("Cannot calibrate against an empty quote set");
    }

    // Standard Nelder-Mead over the 5 Heston parameters. The best vertex is
    // never discarded by any of the four moves, so the recorded best
    // objective is non-increasing by construction - a property the test
    // suite pins.
    constexpr int kDim = 5;
    constexpr int kVerts = kDim + 1;
    constexpr double kAlpha = 1.0;  // Reflection
    constexpr double kGamma = 2.0;  // Expansion
    constexpr double kRho = 0.5;    // Contraction
    constexpr double kSigma = 0.5;  // Shrink

    double simplex[kVerts][kDim];
    double values[kVerts];

    const double guess[kDim] = {initialGuess.v0, initialGuess.theta, initialGuess.kappa,
                                initialGuess.xi, initialGuess.rho};
    for (int v = 0; v < kVerts; ++v) {
        for (int d = 0; d < kDim; ++d) simplex[v][d] = guess[d];
        if (v > 0) {
            const int d = v - 1;
            const double step = (guess[d] != 0.0) ? 0.15 * std::fabs(guess[d]) : 0.05;
            simplex[v][d] = std::clamp(guess[d] + step, kLo[d], kHi[d]);
        }
        values[v] = objective(simplex[v], market);
    }

    Result result{initialGuess, 0.0, 0, {}};
    result.objectiveHistory.reserve(maxIterations_);

    int order[kVerts];
    for (size_t iter = 0; iter < maxIterations_; ++iter) {
        for (int v = 0; v < kVerts; ++v) order[v] = v;
        std::sort(order, order + kVerts,
                  [&](int a, int b) { return values[a] < values[b]; });
        const int best = order[0];
        const int worst = order[kVerts - 1];
        const int secondWorst = order[kVerts - 2];

        result.objectiveHistory.push_back(values[best]);
        if (values[worst] - values[best] < tolerance_) {
            break;
        }

        double centroid[kDim] = {0, 0, 0, 0, 0};
        for (int v = 0; v < kVerts; ++v) {
            if (v == worst) continue;
            for (int d = 0; d < kDim; ++d) centroid[d] += simplex[v][d];
        }
        for (int d = 0; d < kDim; ++d) centroid[d] /= kDim;

        double reflected[kDim];
        for (int d = 0; d < kDim; ++d) {
            reflected[d] = centroid[d] + kAlpha * (centroid[d] - simplex[worst][d]);
        }
        const double fReflected = objective(reflected, market);

        if (fReflected < values[best]) {
            double expanded[kDim];
            for (int d = 0; d < kDim; ++d) {
                expanded[d] = centroid[d] + kGamma * (reflected[d] - centroid[d]);
            }
            const double fExpanded = objective(expanded, market);
            const double* take = (fExpanded < fReflected) ? expanded : reflected;
            const double fTake = std::min(fExpanded, fReflected);
            for (int d = 0; d < kDim; ++d) simplex[worst][d] = take[d];
            values[worst] = fTake;
        } else if (fReflected < values[secondWorst]) {
            for (int d = 0; d < kDim; ++d) simplex[worst][d] = reflected[d];
            values[worst] = fReflected;
        } else {
            double contracted[kDim];
            for (int d = 0; d < kDim; ++d) {
                contracted[d] = centroid[d] + kRho * (simplex[worst][d] - centroid[d]);
            }
            const double fContracted = objective(contracted, market);
            if (fContracted < values[worst]) {
                for (int d = 0; d < kDim; ++d) simplex[worst][d] = contracted[d];
                values[worst] = fContracted;
            } else {
                // Shrink toward the best vertex.
                for (int v = 0; v < kVerts; ++v) {
                    if (v == best) continue;
                    for (int d = 0; d < kDim; ++d) {
                        simplex[v][d] =
                            simplex[best][d] + kSigma * (simplex[v][d] - simplex[best][d]);
                    }
                    values[v] = objective(simplex[v], market);
                }
            }
        }
        ++result.iterations;
    }

    int finalBest = 0;
    for (int v = 1; v < kVerts; ++v) {
        if (values[v] < values[finalBest]) finalBest = v;
    }
    result.params = HestonParams(simplex[finalBest][0], simplex[finalBest][1],
                                 simplex[finalBest][2], simplex[finalBest][3],
                                 simplex[finalBest][4]);
    result.objective = values[finalBest];
    return result;
}

} // namespace Options
