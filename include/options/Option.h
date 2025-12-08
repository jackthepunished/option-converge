#ifndef OPTION_H
#define OPTION_H

#include <string>
#include <stdexcept>

namespace Options {

// Enumeration for option types
enum class OptionType {
    CALL,
    PUT
};

// Enumeration for exercise styles
enum class ExerciseType {
    EUROPEAN,
    AMERICAN
};

// Enumeration for dividend types
enum class DividendType {
    NONE,
    CONTINUOUS,
    DISCRETE
};

// Structure to hold option parameters
struct OptionParams {
    double spotPrice;           // Current stock price (S)
    double strikePrice;         // Strike price (K)
    double riskFreeRate;        // Risk-free interest rate (r)
    double volatility;          // Volatility (sigma)
    double timeToMaturity;      // Time to expiration in years (T)
    double dividendYield;       // Continuous dividend yield (q)
    OptionType optionType;      // Call or Put
    ExerciseType exerciseType;  // European or American

    // Constructor with defaults
    OptionParams(double S = 100.0, double K = 100.0, double r = 0.05,
                 double sigma = 0.20, double T = 1.0, double q = 0.0,
                 OptionType type = OptionType::CALL,
                 ExerciseType exercise = ExerciseType::EUROPEAN)
        : spotPrice(S), strikePrice(K), riskFreeRate(r),
          volatility(sigma), timeToMaturity(T), dividendYield(q),
          optionType(type), exerciseType(exercise) {
        validate();
    }

    // Validation
    void validate() const {
        if (spotPrice <= 0) throw std::invalid_argument("Spot price must be positive");
        if (strikePrice <= 0) throw std::invalid_argument("Strike price must be positive");
        if (volatility < 0) throw std::invalid_argument("Volatility must be non-negative");
        if (timeToMaturity <= 0) throw std::invalid_argument("Time to maturity must be positive");
    }

    // Helper methods
    [[nodiscard]] constexpr bool isCall() const noexcept { return optionType == OptionType::CALL; }
    [[nodiscard]] constexpr bool isPut() const noexcept { return optionType == OptionType::PUT; }
    [[nodiscard]] constexpr bool isEuropean() const noexcept { return exerciseType == ExerciseType::EUROPEAN; }
    [[nodiscard]] constexpr bool isAmerican() const noexcept { return exerciseType == ExerciseType::AMERICAN; }

    [[nodiscard]] std::string toString() const;
};

// Structure to hold Greeks
struct Greeks {
    double delta;    // dV/dS - price sensitivity
    double gamma;    // d²V/dS² - delta sensitivity
    double vega;     // dV/dσ - volatility sensitivity
    double theta;    // dV/dT - time decay
    double rho;      // dV/dr - interest rate sensitivity

    constexpr Greeks() noexcept : delta(0), gamma(0), vega(0), theta(0), rho(0) {}
    constexpr Greeks(double d, double g, double v, double t, double r) noexcept
        : delta(d), gamma(g), vega(v), theta(t), rho(r) {}

    // Rule of Zero - compiler-generated special members are fine
    Greeks(const Greeks&) = default;
    Greeks(Greeks&&) noexcept = default;
    Greeks& operator=(const Greeks&) = default;
    Greeks& operator=(Greeks&&) noexcept = default;
    ~Greeks() = default;

    [[nodiscard]] std::string toString() const;
};

// Structure to hold pricing results
struct PricingResult {
    double price;
    Greeks greeks;
    double standardError;    // For Monte Carlo
    double computationTime;  // In milliseconds
    size_t memoryUsed;       // In bytes

    constexpr PricingResult() noexcept : price(0), greeks(), standardError(0),
                                         computationTime(0), memoryUsed(0) {}

    // Rule of Zero - compiler-generated special members are fine
    PricingResult(const PricingResult&) = default;
    PricingResult(PricingResult&&) noexcept = default;
    PricingResult& operator=(const PricingResult&) = default;
    PricingResult& operator=(PricingResult&&) noexcept = default;
    ~PricingResult() = default;

    [[nodiscard]] std::string toString() const;
};

} // namespace Options

#endif // OPTION_H
