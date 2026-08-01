#ifndef CALIBRATION_H
#define CALIBRATION_H

#include "Heston.h"
#include "Option.h"

#include <cstddef>
#include <vector>

namespace Options {

// One observed market price: a European option identified by strike,
// maturity, and type, with an optional weight for the calibration objective
// (vega weighting or bid-ask weighting in practice; 1.0 by default).
struct MarketQuote {
    double strike;
    double maturity;
    double price;
    OptionType type;
    double weight;

    MarketQuote(double K, double T, double observedPrice, OptionType optType,
                double w = 1.0)
        : strike(K), maturity(T), price(observedPrice), type(optType), weight(w) {
        if (K <= 0.0) throw std::invalid_argument("Quote strike must be positive");
        if (T <= 0.0) throw std::invalid_argument("Quote maturity must be positive");
        if (observedPrice <= 0.0) throw std::invalid_argument("Quote price must be positive");
        if (w <= 0.0) throw std::invalid_argument("Quote weight must be positive");
    }
};

// A quote set plus the observable market state it was taken against.
struct MarketData {
    double spot;
    double riskFreeRate;
    double dividendYield;
    std::vector<MarketQuote> quotes;
};

// Inverts each quote through the Black-Scholes implied volatility solver.
// The resulting (strike, maturity, vol) surface is the market's own
// parameterisation of prices, and the natural target for any smile model.
class ImpliedVolSurface {
public:
    struct Point {
        double strike;
        double maturity;
        double impliedVol;
    };

    [[nodiscard]] static std::vector<Point> build(const MarketData& market);
};

// Fits the five Heston parameters to quoted prices by minimising the
// weighted sum of squared price errors with Nelder-Mead. Derivative-free
// search is the pragmatic choice here: the objective is smooth but its
// gradient in (kappa, xi) is nearly flat, which starves gradient methods,
// and each evaluation is a characteristic-function integration.
class HestonCalibrator {
public:
    struct Result {
        HestonParams params;
        double objective;                      // Final weighted SSE
        size_t iterations;
        std::vector<double> objectiveHistory;  // Best vertex per iteration
    };

    explicit HestonCalibrator(size_t maxIterations = 400, double tolerance = 1e-12)
        : maxIterations_(maxIterations), tolerance_(tolerance) {}

    [[nodiscard]] Result calibrate(const MarketData& market,
                                   const HestonParams& initialGuess) const;

private:
    size_t maxIterations_;
    double tolerance_;

    [[nodiscard]] static double objective(const double p[5], const MarketData& market);
};

} // namespace Options

#endif // CALIBRATION_H
