#include "options/LongstaffSchwartz.h"

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
    if (!params.isAmerican()) {
        throw std::invalid_argument(
            "Longstaff-Schwartz prices American exercise; use MonteCarlo for European options.");
    }

    const size_t N = numPaths_;
    const size_t M = numDates_;
    const double dt = params.timeToMaturity / static_cast<double>(M);
    const double drift =
        (params.riskFreeRate - params.dividendYield - 0.5 * params.volatility * params.volatility) *
        dt;
    const double vol = params.volatility * std::sqrt(dt);
    const double stepDiscount = std::exp(-params.riskFreeRate * dt);

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
                logS += drift + vol * normal(rng);
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
        // Values seen from date k: everything downstream shrinks by one step
        // of discounting, whether or not exercise happens later on the path.
        for (size_t p = 0; p < N; ++p) {
            cashflow[p] *= stepDiscount;
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
        const double value = cashflow[p] * stepDiscount;  // Discount date 1 back to today
        sum += value;
        sumSq += value * value;
    }
    const double n = static_cast<double>(N);
    const double mean = sum / n;
    const double variance = std::fmax(sumSq / n - mean * mean, 0.0);

    // The holder may also exercise immediately, which the date grid
    // (starting at dt) cannot represent; the payoff today floors the price.
    return {std::fmax(mean, intrinsic(params.spotPrice, params)),
            std::sqrt(variance / n)};
}

} // namespace Options
