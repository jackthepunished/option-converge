#include "options/Heston.h"

#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <vector>

namespace Options {

namespace {

using Complex = std::complex<double>;

// Integrand of the P_j probabilities: Re[e^{-iu ln K} phi_j(u) / (iu)], with
// phi_j the "little trap" characteristic function component.
double integrand(double u, int j, const OptionParams& o, const HestonParams& h) {
    const Complex i(0.0, 1.0);
    const double x = std::log(o.spotPrice);
    const double T = o.timeToMaturity;
    const double xi2 = h.xi * h.xi;

    const double uj = (j == 1) ? 0.5 : -0.5;
    const double bj = (j == 1) ? h.kappa - h.rho * h.xi : h.kappa;

    const Complex iu = i * u;
    const Complex beta = bj - h.rho * h.xi * iu;
    const Complex d = std::sqrt(beta * beta - xi2 * (2.0 * uj * iu - u * u));

    // The little trap: work with c = 1/g so the exponent e^{-dT} decays
    // instead of growing, keeping the complex logarithm on its principal
    // branch for all maturities.
    const Complex c = (beta - d) / (beta + d);
    const Complex expMinusDT = std::exp(-d * T);
    const Complex oneMinusCe = 1.0 - c * expMinusDT;

    const Complex D = ((beta - d) / xi2) * ((1.0 - expMinusDT) / oneMinusCe);
    const Complex C = (o.riskFreeRate - o.dividendYield) * iu * T +
                      (h.kappa * h.theta / xi2) *
                          ((beta - d) * T - 2.0 * std::log(oneMinusCe / (1.0 - c)));

    const Complex phi = std::exp(C + D * h.v0 + iu * x);
    return std::real(std::exp(-iu * std::log(o.strikePrice)) * phi / iu);
}

// Adaptive Simpson quadrature. The integrand is smooth and decays like
// 1/u * e^{-const * u^2}, so a modest tolerance resolves it to far better
// than pricing accuracy.
double adaptiveSimpson(double a, double b, double fa, double fm, double fb, double whole,
                       double tol, int depth, int j, const OptionParams& o,
                       const HestonParams& h) {
    const double m = 0.5 * (a + b);
    const double lm = 0.5 * (a + m);
    const double rm = 0.5 * (m + b);
    const double flm = integrand(lm, j, o, h);
    const double frm = integrand(rm, j, o, h);
    const double left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
    const double right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
    if (depth <= 0 || std::fabs(left + right - whole) < 15.0 * tol) {
        return left + right + (left + right - whole) / 15.0;
    }
    return adaptiveSimpson(a, m, fa, flm, fm, left, 0.5 * tol, depth - 1, j, o, h) +
           adaptiveSimpson(m, b, fm, frm, fb, right, 0.5 * tol, depth - 1, j, o, h);
}

double integrate(int j, const OptionParams& o, const HestonParams& h) {
    // The integrand is finite at 0+ and negligible past u ~ 200 for
    // practical parameters; the adaptive refinement handles the shape.
    const double a = 1e-8;
    const double b = 200.0;
    const double fa = integrand(a, j, o, h);
    const double fm = integrand(0.5 * (a + b), j, o, h);
    const double fb = integrand(b, j, o, h);
    const double whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    return adaptiveSimpson(a, b, fa, fm, fb, whole, 1e-9, 40, j, o, h);
}

} // namespace

double HestonModel::probability(const OptionParams& option, int j) const {
    return 0.5 + integrate(j, option, params_) / 3.14159265358979323846;
}

double HestonModel::price(const OptionParams& option) const {
    if (!option.isEuropean()) {
        throw std::invalid_argument("Heston semi-analytic pricing is European only.");
    }
    if (option.hasDiscreteDividends()) {
        throw std::invalid_argument("Heston pricing does not model discrete dividends.");
    }

    const double discS =
        option.spotPrice * std::exp(-option.dividendYield * option.timeToMaturity);
    const double discK =
        option.strikePrice * std::exp(-option.riskFreeRate * option.timeToMaturity);
    const double p1 = probability(option, 1);
    const double p2 = probability(option, 2);

    // Call from the P1/P2 representation; put from the complementary
    // probabilities rather than parity, so parity remains a genuine test of
    // the integration instead of an identity built into the code.
    return option.isCall() ? discS * p1 - discK * p2
                           : discK * (1.0 - p2) - discS * (1.0 - p1);
}

HestonMonteCarlo::HestonMonteCarlo(const HestonParams& params, size_t numPaths,
                                   size_t numSteps, unsigned seed)
    : params_(params), numPaths_(numPaths), numSteps_(numSteps), seed_(seed) {
    if (numPaths_ == 0) throw std::invalid_argument("Heston MC needs at least one path");
    if (numSteps_ == 0) throw std::invalid_argument("Heston MC needs at least one step");
}

HestonMonteCarlo::Result HestonMonteCarlo::price(const OptionParams& option) const {
    if (!option.isEuropean()) {
        throw std::invalid_argument("Heston Monte Carlo prices European options only.");
    }
    if (option.hasDiscreteDividends()) {
        throw std::invalid_argument("Heston Monte Carlo does not model discrete dividends.");
    }

    const double dt = option.timeToMaturity / static_cast<double>(numSteps_);
    const double sqrtDt = std::sqrt(dt);
    const double rho = params_.rho;
    const double rhoBar = std::sqrt(1.0 - rho * rho);
    const double discount = std::exp(-option.riskFreeRate * option.timeToMaturity);
    const double driftBase = (option.riskFreeRate - option.dividendYield) * dt;

    std::mt19937_64 rng(seed_);
    std::normal_distribution<double> normal(0.0, 1.0);

    double sum = 0.0, sumSq = 0.0;
    for (size_t path = 0; path < numPaths_; ++path) {
        double logS = std::log(option.spotPrice);
        double v = params_.v0;
        for (size_t k = 0; k < numSteps_; ++k) {
            // Full truncation: only the positive part of v enters the
            // coefficients; v itself may dip negative and mean-revert back.
            const double vPlus = std::fmax(v, 0.0);
            const double zV = normal(rng);
            const double zS = rho * zV + rhoBar * normal(rng);
            logS += driftBase - 0.5 * vPlus * dt + std::sqrt(vPlus) * sqrtDt * zS;
            v += params_.kappa * (params_.theta - vPlus) * dt +
                 params_.xi * std::sqrt(vPlus) * sqrtDt * zV;
        }
        const double sT = std::exp(logS);
        const double payoff = discount * (option.isCall()
                                              ? std::fmax(sT - option.strikePrice, 0.0)
                                              : std::fmax(option.strikePrice - sT, 0.0));
        sum += payoff;
        sumSq += payoff * payoff;
    }

    const double n = static_cast<double>(numPaths_);
    const double mean = sum / n;
    const double variance = std::fmax(sumSq / n - mean * mean, 0.0);
    return {mean, std::sqrt(variance / n)};
}

namespace {

// Least-squares fit of y on the five-function basis {1, x, x^2, v, v x} via
// the normal equations, Gaussian elimination with partial pivoting. Returns
// false when the system is numerically singular, in which case the caller
// skips the exercise decision at that date.
bool hestonBasisFit(const std::vector<double>& x, const std::vector<double>& v,
                    const std::vector<double>& y, double coeff[5]) {
    constexpr int kDim = 5;
    double xtx[kDim][kDim] = {};
    double xty[kDim] = {};
    const size_t n = x.size();
    for (size_t i = 0; i < n; ++i) {
        const double basis[kDim] = {1.0, x[i], x[i] * x[i], v[i], v[i] * x[i]};
        for (int a = 0; a < kDim; ++a) {
            for (int b = a; b < kDim; ++b) {
                xtx[a][b] += basis[a] * basis[b];
            }
            xty[a] += basis[a] * y[i];
        }
    }
    for (int a = 1; a < kDim; ++a) {
        for (int b = 0; b < a; ++b) {
            xtx[a][b] = xtx[b][a];
        }
    }

    double aug[kDim][kDim + 1];
    for (int a = 0; a < kDim; ++a) {
        for (int b = 0; b < kDim; ++b) aug[a][b] = xtx[a][b];
        aug[a][kDim] = xty[a];
    }
    for (int col = 0; col < kDim; ++col) {
        int pivot = col;
        for (int row = col + 1; row < kDim; ++row) {
            if (std::fabs(aug[row][col]) > std::fabs(aug[pivot][col])) pivot = row;
        }
        if (std::fabs(aug[pivot][col]) < 1e-12) return false;
        if (pivot != col) {
            for (int k = col; k <= kDim; ++k) std::swap(aug[col][k], aug[pivot][k]);
        }
        for (int row = col + 1; row < kDim; ++row) {
            const double factor = aug[row][col] / aug[col][col];
            for (int k = col; k <= kDim; ++k) aug[row][k] -= factor * aug[col][k];
        }
    }
    for (int row = kDim - 1; row >= 0; --row) {
        double sum = aug[row][kDim];
        for (int k = row + 1; k < kDim; ++k) sum -= aug[row][k] * coeff[k];
        coeff[row] = sum / aug[row][row];
    }
    return true;
}

double intrinsicValue(double spot, const OptionParams& o) noexcept {
    return o.isCall() ? std::fmax(spot - o.strikePrice, 0.0)
                      : std::fmax(o.strikePrice - spot, 0.0);
}

} // namespace

HestonLongstaffSchwartz::HestonLongstaffSchwartz(const HestonParams& params, size_t numPaths,
                                                 size_t numSteps, unsigned seed)
    : params_(params), numPaths_(numPaths), numSteps_(numSteps), seed_(seed) {
    if (numPaths_ < 2) throw std::invalid_argument("Heston LSMC needs at least two paths");
    if (numSteps_ == 0) throw std::invalid_argument("Heston LSMC needs at least one step");
}

HestonLongstaffSchwartz::Result HestonLongstaffSchwartz::price(const OptionParams& option) const {
    if (!option.isAmerican()) {
        throw std::invalid_argument(
            "Heston LSMC prices American exercise; use HestonModel or HestonMonteCarlo for "
            "European options.");
    }
    if (option.hasDiscreteDividends()) {
        throw std::invalid_argument("Heston LSMC does not model discrete dividends.");
    }

    const size_t N = numPaths_;
    const size_t M = numSteps_;
    const double dt = option.timeToMaturity / static_cast<double>(M);
    const double sqrtDt = std::sqrt(dt);
    const double rho = params_.rho;
    const double rhoBar = std::sqrt(1.0 - rho * rho);
    const double driftBase = (option.riskFreeRate - option.dividendYield) * dt;
    const double stepDiscount = std::exp(-option.riskFreeRate * dt);

    // Both state variables are needed at every exercise date: the spot for
    // the payoff and the variance for the regression basis. Full truncation
    // Euler, as in the European Heston Monte Carlo.
    std::vector<double> spots(N * M);
    std::vector<double> variances(N * M);
    {
        std::mt19937_64 rng(seed_);
        std::normal_distribution<double> normal(0.0, 1.0);
        for (size_t p = 0; p < N; ++p) {
            double logS = std::log(option.spotPrice);
            double v = params_.v0;
            for (size_t k = 0; k < M; ++k) {
                const double vPlus = std::fmax(v, 0.0);
                const double zV = normal(rng);
                const double zS = rho * zV + rhoBar * normal(rng);
                logS += driftBase - 0.5 * vPlus * dt + std::sqrt(vPlus) * sqrtDt * zS;
                v += params_.kappa * (params_.theta - vPlus) * dt +
                     params_.xi * std::sqrt(vPlus) * sqrtDt * zV;
                spots[p * M + k] = std::exp(logS);
                variances[p * M + k] = std::fmax(v, 0.0);
            }
        }
    }

    std::vector<double> cashflow(N);
    for (size_t p = 0; p < N; ++p) {
        cashflow[p] = intrinsicValue(spots[p * M + (M - 1)], option);
    }

    std::vector<size_t> itmIndex;
    std::vector<double> regX, regV, regY;
    itmIndex.reserve(N);
    regX.reserve(N);
    regV.reserve(N);
    regY.reserve(N);

    for (size_t k = M - 1; k-- > 0;) {
        for (size_t p = 0; p < N; ++p) {
            cashflow[p] *= stepDiscount;
        }

        itmIndex.clear();
        regX.clear();
        regV.clear();
        regY.clear();
        for (size_t p = 0; p < N; ++p) {
            const double exerciseValue = intrinsicValue(spots[p * M + k], option);
            if (exerciseValue > 0.0) {
                itmIndex.push_back(p);
                regX.push_back(spots[p * M + k] / option.strikePrice);
                regV.push_back(variances[p * M + k]);
                regY.push_back(cashflow[p]);
            }
        }
        // Five basis functions need comfortably more support than three.
        if (itmIndex.size() < 25) continue;

        double coeff[5];
        if (!hestonBasisFit(regX, regV, regY, coeff)) continue;

        for (size_t i = 0; i < itmIndex.size(); ++i) {
            const size_t p = itmIndex[i];
            const double x = regX[i];
            const double v = regV[i];
            const double continuation =
                coeff[0] + coeff[1] * x + coeff[2] * x * x + coeff[3] * v + coeff[4] * v * x;
            const double exerciseValue = intrinsicValue(spots[p * M + k], option);
            if (exerciseValue >= continuation) {
                cashflow[p] = exerciseValue;
            }
        }
    }

    double sum = 0.0, sumSq = 0.0;
    for (size_t p = 0; p < N; ++p) {
        const double value = cashflow[p] * stepDiscount;
        sum += value;
        sumSq += value * value;
    }
    const double n = static_cast<double>(N);
    const double mean = sum / n;
    const double variance = std::fmax(sumSq / n - mean * mean, 0.0);

    // Immediate exercise floors the American price, as in the GBM pricer.
    return {std::fmax(mean, intrinsicValue(option.spotPrice, option)),
            std::sqrt(variance / n)};
}

} // namespace Options
