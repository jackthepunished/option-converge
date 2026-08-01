#ifndef HESTON_H
#define HESTON_H

#include "Option.h"

#include <cstddef>

namespace Options {

// Parameters of the Heston (1993) stochastic volatility model:
//   dS = (r - q) S dt + sqrt(v) S dW_S
//   dv = kappa (theta - v) dt + xi sqrt(v) dW_v,   d<W_S, W_v> = rho dt
// The spot's volatility is itself a mean-reverting square-root process, which
// is what produces the smiles and skews constant-volatility Black-Scholes
// cannot.
struct HestonParams {
    double v0;      // Initial variance
    double theta;   // Long-run variance level
    double kappa;   // Mean-reversion speed
    double xi;      // Volatility of variance
    double rho;     // Spot-variance correlation

    HestonParams(double initialVar, double longRunVar, double meanReversion,
                 double volOfVol, double correlation)
        : v0(initialVar), theta(longRunVar), kappa(meanReversion),
          xi(volOfVol), rho(correlation) {
        validate();
    }

    void validate() const {
        if (v0 < 0.0) throw std::invalid_argument("Initial variance must be non-negative");
        if (theta < 0.0) throw std::invalid_argument("Long-run variance must be non-negative");
        if (kappa <= 0.0) throw std::invalid_argument("Mean-reversion speed must be positive");
        if (xi <= 0.0) throw std::invalid_argument("Vol-of-vol must be positive");
        if (rho <= -1.0 || rho >= 1.0) throw std::invalid_argument("Correlation must be in (-1, 1)");
    }

    // The Feller condition keeps the variance process strictly positive.
    // Violating it is legal (and market calibrations routinely do), but the
    // simulation must then handle variance touching zero.
    [[nodiscard]] bool fellerSatisfied() const noexcept {
        return 2.0 * kappa * theta >= xi * xi;
    }
};

// Semi-analytic European pricing by integrating the Heston characteristic
// function, in the Albrecher et al. "little Heston trap" formulation, which
// is numerically stable for long maturities where the original branch-cut
// choice loses continuity. The option's volatility field is ignored; the
// variance dynamics come from HestonParams.
class HestonModel {
public:
    explicit HestonModel(const HestonParams& params) : params_(params) {}

    [[nodiscard]] double price(const OptionParams& option) const;

private:
    HestonParams params_;

    // Risk-neutral exercise probabilities P1 (delta measure) and P2.
    [[nodiscard]] double probability(const OptionParams& option, int j) const;
};

// Monte Carlo pricing under Heston dynamics with the full truncation Euler
// scheme: the variance is floored at zero inside every coefficient but the
// process itself is allowed to go negative and pull back, which is the
// discretisation with the smallest bias among the simple Euler variants.
class HestonMonteCarlo {
public:
    struct Result {
        double price;
        double standardError;
    };

    HestonMonteCarlo(const HestonParams& params, size_t numPaths = 100000,
                     size_t numSteps = 200, unsigned seed = 42);

    [[nodiscard]] Result price(const OptionParams& option) const;

private:
    HestonParams params_;
    size_t numPaths_;
    size_t numSteps_;
    unsigned seed_;
};

} // namespace Options

#endif // HESTON_H
