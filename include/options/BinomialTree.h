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

    // Calculate Greeks using finite differences
    [[nodiscard]] Greeks calculateGreeks(const OptionParams& params) override;

    [[nodiscard]] std::string getName() const noexcept override {
        return "Binomial Tree (CRR, " + std::to_string(steps_) + " steps)";
    }

    // Set number of steps
    void setSteps(size_t steps) noexcept { steps_ = steps; }
    [[nodiscard]] size_t getSteps() const noexcept { return steps_; }

private:
    size_t steps_;  // Number of time steps

    // Calculate option value recursively
    [[nodiscard]] double priceRecursive(const OptionParams& params) const;

    // Calculate option value iteratively (more memory efficient)
    [[nodiscard]] double priceIterative(const OptionParams& params) const;

    // Calculate payoff at expiration
    [[nodiscard]] double payoff(double spotPrice, const OptionParams& params) const noexcept;

    // Check early exercise value for American options
    [[nodiscard]] double earlyExerciseValue(double spotPrice, const OptionParams& params) const noexcept;
};

} // namespace Options

#endif // BINOMIAL_TREE_H
