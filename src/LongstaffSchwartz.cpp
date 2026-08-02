#include "options/LongstaffSchwartz.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace Options {

namespace {

double intrinsic(double spot, const OptionParams& p) noexcept {
    return p.isCall() ? std::fmax(spot - p.strikePrice, 0.0)
                      : std::fmax(p.strikePrice - spot, 0.0);
}

// Least-squares fit of y on {1, x, x^2} via the normal equations, solved with
// Gaussian elimination and partial pivoting. Three unknowns do not justify a
// QR decomposition. Returns false if the system is numerically singular
// (e.g. all regressors nearly identical), in which case the caller skips the
// exercise decision at that date.
bool quadraticFit(const std::vector<double>& x, const std::vector<double>& y,
                  double coeff[3]) {
    double s[5] = {0, 0, 0, 0, 0};  // sums of x^0 .. x^4
    double t[3] = {0, 0, 0};        // sums of y * x^0 .. x^2
    const size_t n = x.size();
    for (size_t i = 0; i < n; ++i) {
        const double xi = x[i];
        const double xi2 = xi * xi;
        s[0] += 1.0;
        s[1] += xi;
        s[2] += xi2;
        s[3] += xi2 * xi;
        s[4] += xi2 * xi2;
        t[0] += y[i];
        t[1] += y[i] * xi;
        t[2] += y[i] * xi2;
    }

    double a[3][4] = {{s[0], s[1], s[2], t[0]},
                      {s[1], s[2], s[3], t[1]},
                      {s[2], s[3], s[4], t[2]}};

    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        for (int row = col + 1; row < 3; ++row) {
            if (std::fabs(a[row][col]) > std::fabs(a[pivot][col])) pivot = row;
        }
        if (std::fabs(a[pivot][col]) < 1e-12) return false;
        if (pivot != col) {
            for (int k = 0; k < 4; ++k) std::swap(a[col][k], a[pivot][k]);
        }
        for (int row = col + 1; row < 3; ++row) {
            const double factor = a[row][col] / a[col][col];
            for (int k = col; k < 4; ++k) a[row][k] -= factor * a[col][k];
        }
    }
    for (int row = 2; row >= 0; --row) {
        double sum = a[row][3];
        for (int k = row + 1; k < 3; ++k) sum -= a[row][k] * coeff[k];
        coeff[row] = sum / a[row][row];
    }
    return true;
}

} // namespace

LongstaffSchwartz::LongstaffSchwartz(size_t numPaths, size_t numExerciseDates, unsigned seed)
    : numPaths_(numPaths), numDates_(numExerciseDates), seed_(seed) {
    if (numPaths_ < 2) throw std::invalid_argument("LSMC needs at least two paths");
    if (numDates_ == 0) throw std::invalid_argument("LSMC needs at least one exercise date");
}

LongstaffSchwartz::Result LongstaffSchwartz::price(const OptionParams& params) const {
    if (params.isEuropean()) {
        throw std::invalid_argument(
            "Longstaff-Schwartz prices American or Bermudan exercise; use MonteCarlo for "
            "European options.");
    }
    if (params.hasDiscreteDividends()) {
        throw std::invalid_argument(
            "Longstaff-Schwartz does not model discrete dividends; use BinomialTree.");
    }

    // The exercise-date grid. American: numDates_ uniform dates ending at
    // expiry. Bermudan: the registered dates, sorted and deduplicated, with
    // expiry appended - the terminal payoff is always available even when
    // no exercise right falls on it.
    std::vector<double> times;
    if (params.isAmerican()) {
        const double dt = params.timeToMaturity / static_cast<double>(numDates_);
        times.reserve(numDates_);
        for (size_t k = 1; k <= numDates_; ++k) {
            times.push_back(dt * static_cast<double>(k));
        }
        times.back() = params.timeToMaturity;
    } else {
        if (params.exerciseDates.empty()) {
            throw std::invalid_argument("Bermudan contract has no exercise dates");
        }
        times = params.exerciseDates;
        std::sort(times.begin(), times.end());
        times.erase(std::unique(times.begin(), times.end(),
                                [](double a, double b) { return std::fabs(a - b) < 1e-12; }),
                    times.end());
        if (params.timeToMaturity - times.back() > 1e-12) {
            times.push_back(params.timeToMaturity);
        } else {
            times.back() = params.timeToMaturity;
        }
    }
    const size_t N = numPaths_;
    const size_t M = times.size();

    // GBM is simulated exactly between consecutive dates, so a sparse
    // Bermudan grid costs nothing in accuracy - only the exercise decisions
    // are discrete, which is the contract's own definition.
    std::vector<double> stepDrift(M), stepVol(M), stepDisc(M);
    {
        const double mu =
            params.riskFreeRate - params.dividendYield - 0.5 * params.volatility * params.volatility;
        double prev = 0.0;
        for (size_t k = 0; k < M; ++k) {
            const double dt = times[k] - prev;
            stepDrift[k] = mu * dt;
            stepVol[k] = params.volatility * std::sqrt(dt);
            stepDisc[k] = std::exp(-params.riskFreeRate * dt);
            prev = times[k];
        }
    }

    // All exercise-date spots are needed for the backward pass, so the full
    // N x M matrix is stored: at 8 bytes a node this is the price of seeing
    // the future that backward induction requires.
    std::vector<double> spots(N * M);
    {
        std::mt19937_64 rng(seed_);
        std::normal_distribution<double> normal(0.0, 1.0);
        for (size_t p = 0; p < N; ++p) {
            double logS = std::log(params.spotPrice);
            for (size_t k = 0; k < M; ++k) {
                logS += stepDrift[k] + stepVol[k] * normal(rng);
                spots[p * M + k] = std::exp(logS);
            }
        }
    }

    // cashflow[p] holds the option value on path p discounted to the date
    // currently being processed; it starts as the terminal payoff.
    std::vector<double> cashflow(N);
    for (size_t p = 0; p < N; ++p) {
        cashflow[p] = intrinsic(spots[p * M + (M - 1)], params);
    }

    std::vector<size_t> itmIndex;
    std::vector<double> regX, regY;
    itmIndex.reserve(N);
    regX.reserve(N);
    regY.reserve(N);

    for (size_t k = M - 1; k-- > 0;) {
        // Values seen from date k: everything downstream shrinks by the
        // discount over (times[k], times[k+1]], whether or not exercise
        // happens later on the path.
        for (size_t p = 0; p < N; ++p) {
            cashflow[p] *= stepDisc[k + 1];
        }

        itmIndex.clear();
        regX.clear();
        regY.clear();
        for (size_t p = 0; p < N; ++p) {
            const double exerciseValue = intrinsic(spots[p * M + k], params);
            if (exerciseValue > 0.0) {
                itmIndex.push_back(p);
                regX.push_back(spots[p * M + k] / params.strikePrice);
                regY.push_back(cashflow[p]);
            }
        }
        // Too few in-the-money paths to support a three-parameter fit.
        if (itmIndex.size() < 10) continue;

        double coeff[3];
        if (!quadraticFit(regX, regY, coeff)) continue;

        for (size_t i = 0; i < itmIndex.size(); ++i) {
            const size_t p = itmIndex[i];
            const double x = regX[i];
            const double continuation = coeff[0] + coeff[1] * x + coeff[2] * x * x;
            const double exerciseValue = intrinsic(spots[p * M + k], params);
            if (exerciseValue >= continuation) {
                cashflow[p] = exerciseValue;
            }
        }
    }

    double sum = 0.0, sumSq = 0.0;
    for (size_t p = 0; p < N; ++p) {
        const double value = cashflow[p] * stepDisc[0];  // Discount the first date back to today
        sum += value;
        sumSq += value * value;
    }
    const double n = static_cast<double>(N);
    const double mean = sum / n;
    const double variance = std::fmax(sumSq / n - mean * mean, 0.0);

    // An American holder may also exercise immediately, which the date grid
    // (starting after 0) cannot represent, so the payoff today floors the
    // price. A Bermudan holder has no such right between dates.
    const double floorValue =
        params.isAmerican() ? intrinsic(params.spotPrice, params) : 0.0;
    return {std::fmax(mean, floorValue), std::sqrt(variance / n)};
}

} // namespace Options
