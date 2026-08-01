#ifndef ASIAN_OPTION_H
#define ASIAN_OPTION_H

#include "Option.h"

#include <cstddef>

namespace Options {

enum class AveragingType {
    ARITHMETIC,
    GEOMETRIC
};

// An Asian contract is a vanilla European payoff applied to the average of
// the spot over discrete fixing dates, uniformly spaced over (0, T]. Wrapped
// around OptionParams like BarrierParams, leaving the engine interface alone.
struct AsianParams {
    OptionParams option;
    AveragingType averaging;
    size_t numFixings;

    AsianParams(const OptionParams& opt, AveragingType avg, size_t fixings)
        : option(opt), averaging(avg), numFixings(fixings) {
        validate();
    }

    void validate() const {
        if (numFixings == 0) throw std::invalid_argument("Asian option needs at least one fixing");
        if (!option.isEuropean()) {
            throw std::invalid_argument("Asian pricing supports European exercise only");
        }
        if (option.hasDiscreteDividends()) {
            // The average samples the path, not just the terminal value, so
            // the escrowed shortcut does not apply.
            throw std::invalid_argument("Asian pricing does not model discrete dividends");
        }
    }
};

// Closed form for the discretely monitored geometric-average Asian option.
// The geometric mean of jointly lognormal fixings is itself lognormal, so
// pricing reduces to the Black-Scholes integral with adjusted drift and
// variance. With a single fixing the contract degenerates to the vanilla
// European option and the formula reproduces Black-Scholes exactly.
class AnalyticGeometricAsian {
public:
    [[nodiscard]] static double price(const AsianParams& params);
};

// Monte Carlo Asian pricer. Simulates the fixings exactly under GBM. For
// arithmetic averaging the geometric payoff on the same path serves as a
// control variate, since its discrete closed form is known exactly and the
// two averages are strongly correlated.
class AsianMonteCarlo {
public:
    struct Result {
        double price;
        double standardError;
    };

    explicit AsianMonteCarlo(size_t numPaths = 100000, unsigned seed = 42);

    [[nodiscard]] Result price(const AsianParams& params, bool useControlVariate = true) const;

private:
    size_t numPaths_;
    unsigned seed_;
};

} // namespace Options

#endif // ASIAN_OPTION_H
