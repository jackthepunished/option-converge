#ifndef BINOMIAL_TREE_H
#define BINOMIAL_TREE_H

#include "PricingEngine.h"
#include <vector>

namespace Options {

// Binomial Tree pricing engine using Cox-Ross-Rubinstein (CRR) parameterization
class BinomialTree final : public PricingEngine {
public:
    // Constructor
    explicit BinomialTree(size_t steps = 100);
    ~BinomialTree() override = default;

    // Allow copying and moving
    BinomialTree(const BinomialTree&) = default;
    BinomialTree(BinomialTree&&) noexcept = default;
    BinomialTree& operator=(const BinomialTree&) = default;
    BinomialTree& operator=(BinomialTree&&) noexcept = default;

    // Price an option using binomial tree
    [[nodiscard]] PricingResult price(const OptionParams& params) override;

    // Delta and Gamma read directly off the tree nodes; Vega, Theta, and Rho
    // by finite differences. Bumping the spot does not work for second-order
    // Greeks on a lattice: the price is piecewise-linear in S, so the default
    // Gamma bump lands inside the noise floor. The step-1 and step-2 node
    // layers, which backward induction produces anyway, ARE the bumped
    // repricings - at exactly the spots the lattice can represent.
    [[nodiscard]] Greeks calculateGreeks(const OptionParams& params) override;

    [[nodiscard]] std::string getName() const noexcept override {
        return "Binomial Tree (CRR, " + std::to_string(steps_) + " steps)";
    }

    // Set number of steps
    void setSteps(size_t steps) noexcept { steps_ = steps; }
    [[nodiscard]] size_t getSteps() const noexcept { return steps_; }

private:
    size_t steps_;  // Number of time steps

    // Node values captured from the early tree layers during backward
    // induction, from which Delta and Gamma follow without repricing.
    struct NodeGreeks {
        double delta;
        double gamma;
        bool valid;  // False when the tree is too shallow (fewer than 2 steps)
    };

    // Calculate option value recursively
    [[nodiscard]] double priceRecursive(const OptionParams& params) const;

    // Calculate option value iteratively (more memory efficient). When
    // nodeGreeks is non-null, Delta and Gamma are extracted on the way.
    [[nodiscard]] double priceIterative(const OptionParams& params,
                                        NodeGreeks* nodeGreeks = nullptr) const;

    // Calculate payoff at expiration
    [[nodiscard]] double payoff(double spotPrice, const OptionParams& params) const noexcept;

    // Check early exercise value for American options
    [[nodiscard]] double earlyExerciseValue(double spotPrice, const OptionParams& params) const noexcept;
};

} // namespace Options

#endif // BINOMIAL_TREE_H
