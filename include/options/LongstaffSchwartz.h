#ifndef LONGSTAFF_SCHWARTZ_H
#define LONGSTAFF_SCHWARTZ_H

#include "Option.h"

#include <cstddef>

namespace Options {

// Longstaff-Schwartz least-squares Monte Carlo for American exercise.
//
// Forward simulation cannot see continuation values, which is why the plain
// MonteCarlo engine refuses American contracts. Longstaff-Schwartz recovers
// them by regression: marching backward over the exercise dates, the
// discounted realised cashflows of in-the-money paths are regressed on a
// polynomial basis in moneyness, and the fitted value stands in for the
// continuation value in the exercise decision. Regression uses only
// in-the-money paths, per the original paper, since only there does the
// exercise decision depend on the fit.
class LongstaffSchwartz {
public:
    struct Result {
        double price;
        double standardError;
    };

    explicit LongstaffSchwartz(size_t numPaths = 100000, size_t numExerciseDates = 50,
                               unsigned seed = 42);

    // Price an American option. Throws std::invalid_argument for European
    // contracts - those belong to the plain MonteCarlo engine.
    [[nodiscard]] Result price(const OptionParams& params) const;

private:
    size_t numPaths_;
    size_t numDates_;
    unsigned seed_;
};

} // namespace Options

#endif // LONGSTAFF_SCHWARTZ_H
