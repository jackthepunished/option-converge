#include "options/BinomialTree.h"

#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Options {

namespace {

// Cox-Ross-Rubinstein lattice parameters for a given set of option params.
struct Lattice {
    double dt;        // Length of one time step, in years
    double u;         // Up factor
    double d;         // Down factor
    double p;         // Risk-neutral probability of an up move
    double discount;  // Per-step discount factor
};

[[nodiscard]] Lattice buildLattice(const OptionParams& params, size_t steps) {
    const double dt = params.timeToMaturity / static_cast<double>(steps);
    const double u = std::exp(params.volatility * std::sqrt(dt));
    const double d = 1.0 / u;
    const double growth = std::exp((params.riskFreeRate - params.dividendYield) * dt);
    const double p = (growth - d) / (u - d);

    // A lattice with p outside [0,1] admits arbitrage and produces meaningless
    // prices. This happens when the step is too coarse for the volatility.
    if (p < 0.0 || p > 1.0) {
        throw std::invalid_argument(
            "Binomial lattice is not arbitrage-free with the given step count; "
            "increase the number of steps.");
    }

    return Lattice{dt, u, d, p, std::exp(-params.riskFreeRate * dt)};
}

} // namespace

BinomialTree::BinomialTree(size_t steps) : steps_(steps) {
    if (steps == 0) {
        throw std::invalid_argument("Binomial tree requires at least one step");
    }
}

// Payoff at a node given the underlying price there.
double BinomialTree::payoff(double spotPrice, const OptionParams& params) const noexcept {
    return params.isCall() ? std::fmax(spotPrice - params.strikePrice, 0.0)
                           : std::fmax(params.strikePrice - spotPrice, 0.0);
}

// Value of exercising immediately. For a vanilla option this is the intrinsic
// value, which is the same expression as the terminal payoff.
double BinomialTree::earlyExerciseValue(double spotPrice, const OptionParams& params) const noexcept {
    return payoff(spotPrice, params);
}

// Backward induction over a single rolling vector of node values.
// Memory is O(steps) rather than the O(steps^2) a full tree would need.
double BinomialTree::priceIterative(const OptionParams& params) const {
    const Lattice lat = buildLattice(params, steps_);
    const size_t n = steps_;

    // Terminal layer: values[i] is the node reached by i up moves.
    std::vector<double> values(n + 1);
    for (size_t i = 0; i <= n; ++i) {
        const double spot = params.spotPrice *
                            std::pow(lat.u, static_cast<double>(i)) *
                            std::pow(lat.d, static_cast<double>(n - i));
        values[i] = payoff(spot, params);
    }

    const bool american = params.isAmerican();

    for (size_t step = n; step-- > 0;) {
        for (size_t i = 0; i <= step; ++i) {
            const double continuation =
                lat.discount * (lat.p * values[i + 1] + (1.0 - lat.p) * values[i]);

            if (american) {
                const double spot = params.spotPrice *
                                    std::pow(lat.u, static_cast<double>(i)) *
                                    std::pow(lat.d, static_cast<double>(step - i));
                values[i] = std::fmax(continuation, earlyExerciseValue(spot, params));
            } else {
                values[i] = continuation;
            }
        }
    }

    return values[0];
}

// Forward recursion with memoization over the recombining lattice. Equivalent
// to priceIterative; kept as a reference implementation and cross-check.
double BinomialTree::priceRecursive(const OptionParams& params) const {
    const Lattice lat = buildLattice(params, steps_);
    const size_t n = steps_;

    // Triangular memo table indexed by (step, upMoves).
    const size_t size = (n + 1) * (n + 2) / 2;
    std::vector<double> memo(size);
    std::vector<char> filled(size, 0);

    const auto index = [](size_t step, size_t up) noexcept {
        return step * (step + 1) / 2 + up;
    };

    const auto spotAt = [&](size_t step, size_t up) noexcept {
        return params.spotPrice * std::pow(lat.u, static_cast<double>(up)) *
               std::pow(lat.d, static_cast<double>(step - up));
    };

    const bool american = params.isAmerican();

    // Explicit stack rather than native recursion: a 10k-step tree would
    // otherwise recurse 10k frames deep.
    struct Frame {
        size_t step;
        size_t up;
        bool expanded;
    };
    std::vector<Frame> stack;
    stack.push_back({0, 0, false});

    while (!stack.empty()) {
        const Frame frame = stack.back();
        stack.pop_back();

        const size_t idx = index(frame.step, frame.up);
        if (filled[idx]) {
            continue;
        }

        if (frame.step == n) {
            memo[idx] = payoff(spotAt(frame.step, frame.up), params);
            filled[idx] = 1;
            continue;
        }

        const size_t downIdx = index(frame.step + 1, frame.up);
        const size_t upIdx = index(frame.step + 1, frame.up + 1);

        if (!frame.expanded && (!filled[downIdx] || !filled[upIdx])) {
            // Revisit this node once both children are resolved.
            stack.push_back({frame.step, frame.up, true});
            if (!filled[upIdx]) {
                stack.push_back({frame.step + 1, frame.up + 1, false});
            }
            if (!filled[downIdx]) {
                stack.push_back({frame.step + 1, frame.up, false});
            }
            continue;
        }

        const double continuation =
            lat.discount * (lat.p * memo[upIdx] + (1.0 - lat.p) * memo[downIdx]);

        memo[idx] = american
                        ? std::fmax(continuation, earlyExerciseValue(spotAt(frame.step, frame.up), params))
                        : continuation;
        filled[idx] = 1;
    }

    return memo[index(0, 0)];
}

PricingResult BinomialTree::price(const OptionParams& params) {
    const auto start = std::chrono::high_resolution_clock::now();

    PricingResult result;
    result.price = priceIterative(params);

    const auto end = std::chrono::high_resolution_clock::now();
    result.computationTime =
        std::chrono::duration<double, std::milli>(end - start).count();

    // Rolling vector of node values is the dominant allocation.
    result.memoryUsed = (steps_ + 1) * sizeof(double);

    return result;
}

Greeks BinomialTree::calculateGreeks(const OptionParams& params) {
    // No closed form on a lattice; fall back to the finite-difference Greeks.
    return PricingEngine::calculateGreeks(params);
}

} // namespace Options
