#include "options/PricingEngine.h"

namespace Options {

// Calculate Greeks using finite differences (default implementation)
Greeks PricingEngine::calculateGreeks(const OptionParams& params) {
    Greeks greeks;

    greeks.delta = calculateDelta(params);
    greeks.gamma = calculateGamma(params);
    greeks.vega = calculateVega(params);
    greeks.theta = calculateTheta(params);
    greeks.rho = calculateRho(params);

    return greeks;
}

// Delta: dV/dS (sensitivity to spot price)
double PricingEngine::calculateDelta(const OptionParams& params, double h) {
    OptionParams paramsUp = params;
    paramsUp.spotPrice += h;
    const double priceUp = price(paramsUp).price;

    OptionParams paramsDown = params;
    paramsDown.spotPrice -= h;
    const double priceDown = price(paramsDown).price;

    return (priceUp - priceDown) / (2.0 * h);
}

// Gamma: d²V/dS² (sensitivity of delta to spot price)
double PricingEngine::calculateGamma(const OptionParams& params, double h) {
    const double priceCenter = price(params).price;

    OptionParams paramsUp = params;
    paramsUp.spotPrice += h;
    const double priceUp = price(paramsUp).price;

    OptionParams paramsDown = params;
    paramsDown.spotPrice -= h;
    const double priceDown = price(paramsDown).price;

    return (priceUp - 2.0 * priceCenter + priceDown) / (h * h);
}

// Vega: dV/dσ (sensitivity to volatility)
double PricingEngine::calculateVega(const OptionParams& params, double h) {
    OptionParams paramsUp = params;
    paramsUp.volatility += h;
    const double priceUp = price(paramsUp).price;

    OptionParams paramsDown = params;
    paramsDown.volatility -= h;
    const double priceDown = price(paramsDown).price;

    // Vega is quoted per 1% change in volatility, matching the analytic engine
    return (priceUp - priceDown) / (2.0 * h) / 100.0;
}

// Theta: dV/dT (time decay)
double PricingEngine::calculateTheta(const OptionParams& params, double h) {
    // Ensure we don't go negative
    if (params.timeToMaturity <= h) {
        return 0.0;
    }

    const double priceNow = price(params).price;

    OptionParams paramsForward = params;
    paramsForward.timeToMaturity -= h;
    const double priceForward = price(paramsForward).price;

    // Theta is negative for long positions (time decay)
    // Return daily theta (per day), matching the analytic engine
    return (priceForward - priceNow) / h / 365.0;
}

// Rho: dV/dr (sensitivity to interest rate)
double PricingEngine::calculateRho(const OptionParams& params, double h) {
    OptionParams paramsUp = params;
    paramsUp.riskFreeRate += h;
    const double priceUp = price(paramsUp).price;

    OptionParams paramsDown = params;
    paramsDown.riskFreeRate -= h;
    const double priceDown = price(paramsDown).price;

    // Rho is quoted per 1% change in interest rate, matching the analytic engine
    return (priceUp - priceDown) / (2.0 * h) / 100.0;
}

} // namespace Options
