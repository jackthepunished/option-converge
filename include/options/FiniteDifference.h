#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include "PricingEngine.h"

namespace Options {

// Discretization schemes for the Black-Scholes PDE, expressed as the theta
// parameter of a unified theta-scheme: explicit (theta = 0, conditionally
// stable), fully implicit (theta = 1, unconditionally stable, O(dt)),
// Crank-Nicolson (theta = 1/2, unconditionally stable, O(dt^2)).
enum class FDScheme {
    EXPLICIT,
    IMPLICIT,
    CRANK_NICOLSON
};

// Finite difference pricing engine. Solves the Black-Scholes PDE backward in
// time on a uniform log-spot grid, which makes the PDE coefficients constant
// so a single tridiagonal stencil serves every time step. Supports European
// and American exercise; American values apply the early-exercise projection
// max(continuation, intrinsic) after each step.
class FiniteDifference final : public PricingEngine {
public:
    explicit FiniteDifference(size_t spotSteps = 200, size_t timeSteps = 200,
                              FDScheme scheme = FDScheme::CRANK_NICOLSON);
    ~FiniteDifference() override = default;

    FiniteDifference(const FiniteDifference&) = default;
    FiniteDifference(FiniteDifference&&) noexcept = default;
    FiniteDifference& operator=(const FiniteDifference&) = default;
    FiniteDifference& operator=(FiniteDifference&&) noexcept = default;

    [[nodiscard]] PricingResult price(const OptionParams& params) override;

    [[nodiscard]] std::string getName() const override;

    void setScheme(FDScheme scheme) noexcept { scheme_ = scheme; }
    [[nodiscard]] FDScheme getScheme() const noexcept { return scheme_; }
    [[nodiscard]] size_t getSpotSteps() const noexcept { return spotSteps_; }
    [[nodiscard]] size_t getTimeSteps() const noexcept { return timeSteps_; }

private:
    size_t spotSteps_;   // Grid intervals in log-spot (kept even so the spot lands on a node)
    size_t timeSteps_;   // Grid intervals in time (raised automatically for explicit stability)
    FDScheme scheme_;

    [[nodiscard]] double solveGrid(const OptionParams& params, size_t& memoryUsed) const;

    [[nodiscard]] static double payoff(double spot, const OptionParams& params) noexcept;

    // Dirichlet boundary values at the grid edges for time-to-maturity tau.
    [[nodiscard]] static double lowerBoundary(double spot, double tau,
                                              const OptionParams& params) noexcept;
    [[nodiscard]] static double upperBoundary(double spot, double tau,
                                              const OptionParams& params) noexcept;
};

} // namespace Options

#endif // FINITE_DIFFERENCE_H
