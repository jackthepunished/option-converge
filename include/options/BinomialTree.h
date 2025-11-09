#ifndef BINOMIAL_TREE_H
#define BINOMIAL_TREE_H

#include "PricingEngine.h"
#include <vector>

namespace Options {

// Binomial Tree pricing engine using Cox-Ross-Rubinstein (CRR) parameterization
class BinomialTree : public PricingEngine {
public:
    // Constructor
    explicit BinomialTree(size_t steps = 100);

    // Price an option using binomial tree
    PricingResult price(const OptionParams& params) override;

    // Calculate Greeks using finite differences
    Greeks calculateGreeks(const OptionParams& params) override;

    std::string getName() const override {
        return "Binomial Tree (CRR, " + std::to_string(steps_) + " steps)";
    }

    // Set number of steps
    void setSteps(size_t steps) { steps_ = steps; }
    size_t getSteps() const { return steps_; }

private:
    size_t steps_;  // Number of time steps

    // Calculate option value recursively
    double priceRecursive(const OptionParams& params) const;

    // Calculate option value iteratively (more memory efficient)
    double priceIterative(const OptionParams& params) const;

    // Calculate payoff at expiration
    double payoff(double spotPrice, const OptionParams& params) const;

    // Check early exercise value for American options
    double earlyExerciseValue(double spotPrice, const OptionParams& params) const;
};

} // namespace Options

#endif // BINOMIAL_TREE_H
