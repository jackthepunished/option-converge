#ifndef DIGITAL_OPTION_H
#define DIGITAL_OPTION_H

#include "Option.h"

#include <cstddef>

namespace Options {

// Cash-or-nothing pays a fixed amount if the option finishes in the money;
// asset-or-nothing pays the asset itself. Vanillas decompose into them:
//   call = asset-or-nothing call - K * cash-or-nothing call (payout 1)
// which the test suite pins as an exact identity.
enum class DigitalType {
    CASH_OR_NOTHING,
    ASSET_OR_NOTHING
};

struct DigitalParams {
    OptionParams option;
    DigitalType digitalType;
    double cashPayout;  // Paid by CASH_OR_NOTHING; ignored by ASSET_OR_NOTHING

    DigitalParams(const OptionParams& opt, DigitalType type, double payout = 1.0)
        : option(opt), digitalType(type), cashPayout(payout) {
        validate();
    }

    void validate() const {
        if (!option.isEuropean()) {
            throw std::invalid_argument("Digital pricing supports European exercise only");
        }
        if (option.hasDiscreteDividends()) {
            throw std::invalid_argument("Digital pricing does not model discrete dividends");
        }
        if (digitalType == DigitalType::CASH_OR_NOTHING && cashPayout <= 0.0) {
            throw std::invalid_argument("Cash payout must be positive");
        }
    }
};

// Closed forms under Black-Scholes. The digitals are the two halves of the
// vanilla formula priced separately: Q e^{-rT} N(d2) and S e^{-qT} N(d1).
class AnalyticDigital {
public:
    [[nodiscard]] static double price(const DigitalParams& params);

    // Analytic delta. Digitals concentrate all their optionality at the
    // strike: near expiry the delta at the money grows like 1/(sigma
    // sqrt(T)), the spike that makes them notoriously hard to hedge.
    [[nodiscard]] static double delta(const DigitalParams& params);
};

// CRR lattice pricing of the indicator payoff. Deliberately included for
// what it teaches: a discontinuous payoff destroys the lattice's smooth
// O(1/N) convergence - the error oscillates as the strike moves between
// terminal nodes and decays only like O(1/sqrt(N)). The README documents
// this; the tests bound it loosely rather than pretending it away.
class DigitalLattice {
public:
    explicit DigitalLattice(size_t steps = 1000);

    [[nodiscard]] double price(const DigitalParams& params) const;

private:
    size_t steps_;
};

} // namespace Options

#endif // DIGITAL_OPTION_H
