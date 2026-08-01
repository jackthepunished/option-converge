#include "options/Heston.h"

#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>

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

} // namespace Options
