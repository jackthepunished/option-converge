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
    OptionParams paramsDown = params;

    paramsUp.spotPrice = params.spotPrice + h;
    paramsDown.spotPrice = params.spotPrice - h;

    double priceUp = price(paramsUp).price;
    double priceDown = price(paramsDown).price;

    return (priceUp - priceDown) / (2.0 * h);
}

// Gamma: d²V/dS² (sensitivity of delta to spot price)
double PricingEngine::calculateGamma(const OptionParams& params, double h) {
    OptionParams paramsUp = params;
    OptionParams paramsDown = params;
    OptionParams paramsCenter = params;

    paramsUp.spotPrice = params.spotPrice + h;
    paramsDown.spotPrice = params.spotPrice - h;

    double priceUp = price(paramsUp).price;
    double priceDown = price(paramsDown).price;
    double priceCenter = price(paramsCenter).price;

    return (priceUp - 2.0 * priceCenter + priceDown) / (h * h);
}

// Vega: dV/dσ (sensitivity to volatility)
double PricingEngine::calculateVega(const OptionParams& params, double h) {
    OptionParams paramsUp = params;
    OptionParams paramsDown = params;

    paramsUp.volatility = params.volatility + h;
    paramsDown.volatility = params.volatility - h;

    double priceUp = price(paramsUp).price;
    double priceDown = price(paramsDown).price;

    // Vega is typically quoted per 1% change in volatility
    return (priceUp - priceDown) / (2.0 * h);
}

// Theta: dV/dT (time decay)
double PricingEngine::calculateTheta(const OptionParams& params, double h) {
    OptionParams paramsForward = params;

    // Move time forward by h (reduce time to maturity)
    paramsForward.timeToMaturity = params.timeToMaturity - h;

    // Ensure we don't go negative
    if (paramsForward.timeToMaturity <= 0) {
        return 0.0;
    }

    double priceNow = price(params).price;
    double priceForward = price(paramsForward).price;

    // Theta is negative for long positions (time decay)
    // Return daily theta (per day)
    return (priceForward - priceNow) / h;
}

// Rho: dV/dr (sensitivity to interest rate)
double PricingEngine::calculateRho(const OptionParams& params, double h) {
    OptionParams paramsUp = params;
    OptionParams paramsDown = params;

    paramsUp.riskFreeRate = params.riskFreeRate + h;
    paramsDown.riskFreeRate = params.riskFreeRate - h;

    double priceUp = price(paramsUp).price;
    double priceDown = price(paramsDown).price;

    // Rho is typically quoted per 1% change in interest rate
    return (priceUp - priceDown) / (2.0 * h);
}

} // namespace Options
