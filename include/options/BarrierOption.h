#ifndef BARRIER_OPTION_H
#define BARRIER_OPTION_H

#include "Option.h"

#include <cstddef>

namespace Options {

// The four single-barrier variants. "Out" options die when the barrier is
// touched; "in" options only come alive then. With zero rebate the two are
// complementary: knock-in + knock-out = vanilla, path by path.
enum class BarrierType {
    UP_AND_IN,
    UP_AND_OUT,
    DOWN_AND_IN,
    DOWN_AND_OUT
};

// A barrier contract is a vanilla European option plus a barrier level and
// direction. Kept as a wrapper around OptionParams rather than new fields on
// it, so the vanilla engines' interface is untouched.
struct BarrierParams {
    OptionParams option;
    BarrierType barrierType;
    double barrier;

    BarrierParams(const OptionParams& opt, BarrierType type, double barrierLevel)
        : option(opt), barrierType(type), barrier(barrierLevel) {
        validate();
    }

    void validate() const {
        if (barrier <= 0.0) throw std::invalid_argument("Barrier level must be positive");
        if (!option.isEuropean()) {
            throw std::invalid_argument("Barrier pricing supports European exercise only");
        }
    }

    [[nodiscard]] constexpr bool isUp() const noexcept {
        return barrierType == BarrierType::UP_AND_IN || barrierType == BarrierType::UP_AND_OUT;
    }
    [[nodiscard]] constexpr bool isKnockIn() const noexcept {
        return barrierType == BarrierType::UP_AND_IN || barrierType == BarrierType::DOWN_AND_IN;
    }

    // True if the barrier is already touched at inception, in which case a
    // knock-out is worthless and a knock-in is simply the vanilla option.
    [[nodiscard]] bool alreadyBreached() const noexcept {
        return isUp() ? option.spotPrice >= barrier : option.spotPrice <= barrier;
    }
};

// Closed-form prices for continuously monitored European barrier options
// under Black-Scholes dynamics (Reiner-Rubinstein / Merton). Knock-ins are
// evaluated directly; knock-outs follow from in-out parity against the
// vanilla price, which holds exactly at zero rebate.
class AnalyticBarrier {
public:
    [[nodiscard]] static double price(const BarrierParams& params);
};

// Monte Carlo barrier pricer. Simulates GBM paths exactly at each step and
// corrects for what happens between steps with the Brownian-bridge crossing
// probability, so coarse discrete monitoring approximates the continuously
// monitored contract the analytic formulas price.
class BarrierMonteCarlo {
public:
    struct Result {
        double price;
        double standardError;
    };

    explicit BarrierMonteCarlo(size_t numPaths = 100000, size_t numSteps = 64,
                               unsigned seed = 42);

    [[nodiscard]] Result price(const BarrierParams& params) const;

private:
    size_t numPaths_;
    size_t numSteps_;
    unsigned seed_;
};

} // namespace Options

#endif // BARRIER_OPTION_H
