#include "options/BlackScholes.h"
#include <cmath>
#include <chrono>
#include <stdexcept>

namespace Options {

// Price an option using Black-Scholes formula
PricingResult BlackScholes::price(const OptionParams& params) {
    auto start = std::chrono::high_resolution_clock::now();

    // Validate that this is a European option
    if (!params.isEuropean()) {
        throw std::invalid_argument(
            "Black-Scholes only prices European options. Use Binomial Tree for American options.");
    }

    PricingResult result;

    // Calculate price
    if (params.isCall()) {
        result.price = callPrice(params);
    } else {
        result.price = putPrice(params);
    }

    // Calculate Greeks analytically
    result.greeks = calculateGreeks(params);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    result.computationTime = duration.count() / 1000.0; // Convert to milliseconds

    return result;
}

// Calculate Greeks analytically (overrides base class finite difference method)
Greeks BlackScholes::calculateGreeks(const OptionParams& params) {
    Greeks greeks;

    if (params.isCall()) {
        greeks.delta = callDelta(params);
        greeks.theta = callTheta(params);
        greeks.rho = callRho(params);
    } else {
        greeks.delta = putDelta(params);
        greeks.theta = putTheta(params);
        greeks.rho = putRho(params);
    }

    // Gamma and Vega are the same for calls and puts
    greeks.gamma = gamma(params);
    greeks.vega = vega(params);

    return greeks;
}

// Calculate d1 parameter for Black-Scholes
double BlackScholes::d1(const OptionParams& params) const {
    double S = params.spotPrice;
    double K = params.strikePrice;
    double r = params.riskFreeRate;
    double q = params.dividendYield;
    double sigma = params.volatility;
    double T = params.timeToMaturity;

    return (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
}

// Calculate d2 parameter for Black-Scholes
double BlackScholes::d2(const OptionParams& params) const {
    double sigma = params.volatility;
    double T = params.timeToMaturity;

    return d1(params) - sigma * std::sqrt(T);
}

// Standard normal cumulative distribution function
double BlackScholes::normalCDF(double x) const {
    // Using approximation for standard normal CDF
    // More accurate implementation would use erf() function
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// Standard normal probability density function
double BlackScholes::normalPDF(double x) const {
    return (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * x * x);
}

// Calculate call option price
double BlackScholes::callPrice(const OptionParams& params) const {
    double S = params.spotPrice;
    double K = params.strikePrice;
    double r = params.riskFreeRate;
    double q = params.dividendYield;
    double T = params.timeToMaturity;

    double d1_val = d1(params);
    double d2_val = d2(params);

    double Nd1 = normalCDF(d1_val);
    double Nd2 = normalCDF(d2_val);

    return S * std::exp(-q * T) * Nd1 - K * std::exp(-r * T) * Nd2;
}

// Calculate put option price
double BlackScholes::putPrice(const OptionParams& params) const {
    double S = params.spotPrice;
    double K = params.strikePrice;
    double r = params.riskFreeRate;
    double q = params.dividendYield;
    double T = params.timeToMaturity;

    double d1_val = d1(params);
    double d2_val = d2(params);

    double Nminusd1 = normalCDF(-d1_val);
    double Nminusd2 = normalCDF(-d2_val);

    return K * std::exp(-r * T) * Nminusd2 - S * std::exp(-q * T) * Nminusd1;
}

// Calculate call delta
double BlackScholes::callDelta(const OptionParams& params) const {
    double q = params.dividendYield;
    double T = params.timeToMaturity;
    double d1_val = d1(params);

    return std::exp(-q * T) * normalCDF(d1_val);
}

// Calculate put delta
double BlackScholes::putDelta(const OptionParams& params) const {
    double q = params.dividendYield;
    double T = params.timeToMaturity;
    double d1_val = d1(params);

    return -std::exp(-q * T) * normalCDF(-d1_val);
}

// Calculate gamma (same for calls and puts)
double BlackScholes::gamma(const OptionParams& params) const {
    double S = params.spotPrice;
    double q = params.dividendYield;
    double sigma = params.volatility;
    double T = params.timeToMaturity;
    double d1_val = d1(params);

    return std::exp(-q * T) * normalPDF(d1_val) / (S * sigma * std::sqrt(T));
}

// Calculate vega (same for calls and puts)
double BlackScholes::vega(const OptionParams& params) const {
    double S = params.spotPrice;
    double q = params.dividendYield;
    double T = params.timeToMaturity;
    double d1_val = d1(params);

    // Vega is typically quoted per 1% change in volatility
    // So we divide by 100 to get the sensitivity per 0.01 change
    return S * std::exp(-q * T) * normalPDF(d1_val) * std::sqrt(T) / 100.0;
}

// Calculate call theta
double BlackScholes::callTheta(const OptionParams& params) const {
    double S = params.spotPrice;
    double K = params.strikePrice;
    double r = params.riskFreeRate;
    double q = params.dividendYield;
    double sigma = params.volatility;
    double T = params.timeToMaturity;

    double d1_val = d1(params);
    double d2_val = d2(params);

    double term1 = -(S * normalPDF(d1_val) * sigma * std::exp(-q * T)) / (2.0 * std::sqrt(T));
    double term2 = r * K * std::exp(-r * T) * normalCDF(d2_val);
    double term3 = -q * S * std::exp(-q * T) * normalCDF(d1_val);

    // Return daily theta (divide by 365)
    return (term1 + term2 + term3) / 365.0;
}

// Calculate put theta
double BlackScholes::putTheta(const OptionParams& params) const {
    double S = params.spotPrice;
    double K = params.strikePrice;
    double r = params.riskFreeRate;
    double q = params.dividendYield;
    double sigma = params.volatility;
    double T = params.timeToMaturity;

    double d1_val = d1(params);
    double d2_val = d2(params);

    double term1 = -(S * normalPDF(d1_val) * sigma * std::exp(-q * T)) / (2.0 * std::sqrt(T));
    double term2 = -r * K * std::exp(-r * T) * normalCDF(-d2_val);
    double term3 = q * S * std::exp(-q * T) * normalCDF(-d1_val);

    // Return daily theta (divide by 365)
    return (term1 + term2 + term3) / 365.0;
}

// Calculate call rho
double BlackScholes::callRho(const OptionParams& params) const {
    double K = params.strikePrice;
    double r = params.riskFreeRate;
    double T = params.timeToMaturity;
    double d2_val = d2(params);

    // Rho is typically quoted per 1% change in interest rate
    // So we divide by 100 to get the sensitivity per 0.01 change
    return K * T * std::exp(-r * T) * normalCDF(d2_val) / 100.0;
}

// Calculate put rho
double BlackScholes::putRho(const OptionParams& params) const {
    double K = params.strikePrice;
    double r = params.riskFreeRate;
    double T = params.timeToMaturity;
    double d2_val = d2(params);

    // Rho is typically quoted per 1% change in interest rate
    // So we divide by 100 to get the sensitivity per 0.01 change
    return -K * T * std::exp(-r * T) * normalCDF(-d2_val) / 100.0;
}

} // namespace Options
