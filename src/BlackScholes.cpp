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

    // Pre-calculate all values once for both price and Greeks
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;
    const double sqrtT = std::sqrt(T);

    const double d1_val = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2_val = d1_val - sigma * sqrtT;

    const double Nd1 = normalCDF(d1_val);
    const double Nd2 = normalCDF(d2_val);
    const double nd1 = normalPDF(d1_val);
    const double expMinusQT = std::exp(-q * T);
    const double expMinusRT = std::exp(-r * T);

    // Calculate price using cached values
    if (params.isCall()) {
        result.price = S * expMinusQT * Nd1 - K * expMinusRT * Nd2;

        // Calculate Greeks for call
        result.greeks.delta = expMinusQT * Nd1;
        result.greeks.theta = (-(S * nd1 * sigma * expMinusQT) / (2.0 * sqrtT)
                              + r * K * expMinusRT * Nd2
                              - q * S * expMinusQT * Nd1) / 365.0;
        result.greeks.rho = K * T * expMinusRT * Nd2 / 100.0;
    } else {
        const double NminusD1 = normalCDF(-d1_val);
        const double NminusD2 = normalCDF(-d2_val);

        result.price = K * expMinusRT * NminusD2 - S * expMinusQT * NminusD1;

        // Calculate Greeks for put
        result.greeks.delta = -expMinusQT * NminusD1;
        result.greeks.theta = (-(S * nd1 * sigma * expMinusQT) / (2.0 * sqrtT)
                              - r * K * expMinusRT * NminusD2
                              + q * S * expMinusQT * NminusD1) / 365.0;
        result.greeks.rho = -K * T * expMinusRT * NminusD2 / 100.0;
    }

    // Gamma and Vega are the same for calls and puts
    result.greeks.gamma = expMinusQT * nd1 / (S * sigma * sqrtT);
    result.greeks.vega = S * expMinusQT * nd1 * sqrtT / 100.0;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    result.computationTime = duration.count() / 1000.0; // Convert to milliseconds

    return result;
}

// Calculate Greeks analytically (overrides base class finite difference method)
Greeks BlackScholes::calculateGreeks(const OptionParams& params) {
    Greeks greeks;

    // Pre-calculate common values once
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;
    const double sqrtT = std::sqrt(T);

    const double d1_val = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2_val = d1_val - sigma * sqrtT;

    const double Nd1 = normalCDF(d1_val);
    const double Nd2 = normalCDF(d2_val);
    const double nd1 = normalPDF(d1_val);
    const double expMinusQT = std::exp(-q * T);
    const double expMinusRT = std::exp(-r * T);

    // Calculate all Greeks using cached values
    if (params.isCall()) {
        greeks.delta = expMinusQT * Nd1;
        greeks.theta = (-(S * nd1 * sigma * expMinusQT) / (2.0 * sqrtT)
                       + r * K * expMinusRT * Nd2
                       - q * S * expMinusQT * Nd1) / 365.0;
        greeks.rho = K * T * expMinusRT * Nd2 / 100.0;
    } else {
        greeks.delta = -expMinusQT * normalCDF(-d1_val);
        greeks.theta = (-(S * nd1 * sigma * expMinusQT) / (2.0 * sqrtT)
                       - r * K * expMinusRT * normalCDF(-d2_val)
                       + q * S * expMinusQT * normalCDF(-d1_val)) / 365.0;
        greeks.rho = -K * T * expMinusRT * normalCDF(-d2_val) / 100.0;
    }

    // Gamma and Vega are the same for calls and puts
    greeks.gamma = expMinusQT * nd1 / (S * sigma * sqrtT);
    greeks.vega = S * expMinusQT * nd1 * sqrtT / 100.0;

    return greeks;
}

// Calculate d1 parameter for Black-Scholes
inline double BlackScholes::d1(const OptionParams& params) const noexcept {
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;

    return (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
}

// Calculate d2 parameter for Black-Scholes
inline double BlackScholes::d2(const OptionParams& params) const noexcept {
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;

    return d1(params) - sigma * std::sqrt(T);
}

// Standard normal cumulative distribution function
inline double BlackScholes::normalCDF(double x) const noexcept {
    // Using approximation for standard normal CDF
    // More accurate implementation would use erf() function
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// Standard normal probability density function
inline double BlackScholes::normalPDF(double x) const noexcept {
    constexpr double inv_sqrt_2pi = 0.3989422804014327; // 1/sqrt(2π) precomputed
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

// Calculate call option price
double BlackScholes::callPrice(const OptionParams& params) const {
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double T = params.timeToMaturity;

    const double d1_val = d1(params);
    const double d2_val = d2(params);

    const double Nd1 = normalCDF(d1_val);
    const double Nd2 = normalCDF(d2_val);

    return S * std::exp(-q * T) * Nd1 - K * std::exp(-r * T) * Nd2;
}

// Calculate put option price
double BlackScholes::putPrice(const OptionParams& params) const {
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double T = params.timeToMaturity;

    const double d1_val = d1(params);
    const double d2_val = d2(params);

    const double Nminusd1 = normalCDF(-d1_val);
    const double Nminusd2 = normalCDF(-d2_val);

    return K * std::exp(-r * T) * Nminusd2 - S * std::exp(-q * T) * Nminusd1;
}

// Calculate call delta
double BlackScholes::callDelta(const OptionParams& params) const {
    const double q = params.dividendYield;
    const double T = params.timeToMaturity;
    const double d1_val = d1(params);

    return std::exp(-q * T) * normalCDF(d1_val);
}

// Calculate put delta
double BlackScholes::putDelta(const OptionParams& params) const {
    const double q = params.dividendYield;
    const double T = params.timeToMaturity;
    const double d1_val = d1(params);

    return -std::exp(-q * T) * normalCDF(-d1_val);
}

// Calculate gamma (same for calls and puts)
double BlackScholes::gamma(const OptionParams& params) const {
    const double S = params.spotPrice;
    const double q = params.dividendYield;
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;
    const double d1_val = d1(params);

    return std::exp(-q * T) * normalPDF(d1_val) / (S * sigma * std::sqrt(T));
}

// Calculate vega (same for calls and puts)
double BlackScholes::vega(const OptionParams& params) const {
    const double S = params.spotPrice;
    const double q = params.dividendYield;
    const double T = params.timeToMaturity;
    const double d1_val = d1(params);

    // Vega is typically quoted per 1% change in volatility
    // So we divide by 100 to get the sensitivity per 0.01 change
    return S * std::exp(-q * T) * normalPDF(d1_val) * std::sqrt(T) / 100.0;
}

// Calculate call theta
double BlackScholes::callTheta(const OptionParams& params) const {
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;

    const double d1_val = d1(params);
    const double d2_val = d2(params);

    const double term1 = -(S * normalPDF(d1_val) * sigma * std::exp(-q * T)) / (2.0 * std::sqrt(T));
    const double term2 = r * K * std::exp(-r * T) * normalCDF(d2_val);
    const double term3 = -q * S * std::exp(-q * T) * normalCDF(d1_val);

    // Return daily theta (divide by 365)
    return (term1 + term2 + term3) / 365.0;
}

// Calculate put theta
double BlackScholes::putTheta(const OptionParams& params) const {
    const double S = params.spotPrice;
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double q = params.dividendYield;
    const double sigma = params.volatility;
    const double T = params.timeToMaturity;

    const double d1_val = d1(params);
    const double d2_val = d2(params);

    const double term1 = -(S * normalPDF(d1_val) * sigma * std::exp(-q * T)) / (2.0 * std::sqrt(T));
    const double term2 = -r * K * std::exp(-r * T) * normalCDF(-d2_val);
    const double term3 = q * S * std::exp(-q * T) * normalCDF(-d1_val);

    // Return daily theta (divide by 365)
    return (term1 + term2 + term3) / 365.0;
}

// Calculate call rho
double BlackScholes::callRho(const OptionParams& params) const {
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double T = params.timeToMaturity;
    const double d2_val = d2(params);

    // Rho is typically quoted per 1% change in interest rate
    // So we divide by 100 to get the sensitivity per 0.01 change
    return K * T * std::exp(-r * T) * normalCDF(d2_val) / 100.0;
}

// Calculate put rho
double BlackScholes::putRho(const OptionParams& params) const {
    const double K = params.strikePrice;
    const double r = params.riskFreeRate;
    const double T = params.timeToMaturity;
    const double d2_val = d2(params);

    // Rho is typically quoted per 1% change in interest rate
    // So we divide by 100 to get the sensitivity per 0.01 change
    return -K * T * std::exp(-r * T) * normalCDF(-d2_val) / 100.0;
}

} // namespace Options
