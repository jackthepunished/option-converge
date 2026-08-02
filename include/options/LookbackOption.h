#ifndef LOOKBACK_OPTION_H
#define LOOKBACK_OPTION_H

#include "Option.h"

#include <cstddef>

namespace Options {

// Floating strike: the strike IS the realised extreme (min for calls, max
// for puts), so the option always finishes in the money. Fixed strike: a
// vanilla strike applied to the realised extreme instead of the terminal
// price. Both are priced from inception, with the running extreme starting
// at the current spot.
enum class LookbackType {
    FLOATING_STRIKE,
    FIXED_STRIKE
};

struct LookbackParams {
    OptionParams option;  // strikePrice is ignored for FLOATING_STRIKE
    LookbackType lookbackType;

    LookbackParams(const OptionParams& opt, LookbackType type)
        : option(opt), lookbackType(type) {
        validate();
    }

    void validate() const {
        if (!option.isEuropean()) {
            throw std::invalid_argument("Lookback pricing supports European exercise only");
        }
        if (option.hasDiscreteDividends()) {
            // The extreme samples the whole path, which the escrowed
            // adjustment distorts.
            throw std::invalid_argument("Lookback pricing does not model discrete dividends");
        }
    }
};

// Closed forms for continuously monitored lookbacks under Black-Scholes:
// Goldman-Sosin-Gatto for floating strikes, Conze-Viswanathan for fixed.
// The formulas carry a sigma^2 / (2(r - q)) factor, so the r == q limit is
// a removable singularity this implementation does not take; it refuses
// rather than returning the wrong branch.
class AnalyticLookback {
public:
    [[nodiscard]] static double price(const LookbackParams& params);
};

// Monte Carlo lookback pricer. GBM nodes are simulated exactly and the
// extreme BETWEEN nodes is sampled from the Brownian bridge's known extreme
// distribution, so a coarse time grid prices the continuously monitored
// contract instead of the biased discretely monitored one - the same idea
// as the barrier pricer's crossing correction, applied to the extreme value
// itself.
class LookbackMonteCarlo {
public:
    struct Result {
        double price;
        double standardError;
    };

    explicit LookbackMonteCarlo(size_t numPaths = 100000, size_t numSteps = 64,
                                unsigned seed = 42);

    [[nodiscard]] Result price(const LookbackParams& params) const;

private:
    size_t numPaths_;
    size_t numSteps_;
    unsigned seed_;
};

} // namespace Options

#endif // LOOKBACK_OPTION_H
