#ifndef OPTION_H
#define OPTION_H

#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>

namespace Options {

// Enumeration for option types
enum class OptionType {
    CALL,
    PUT
};

// Enumeration for exercise styles. Bermudan sits between the other two:
// exercise is allowed only on a discrete set of dates, registered via
// OptionParams::addExerciseDate.
enum class ExerciseType {
    EUROPEAN,
    AMERICAN,
    BERMUDAN
};

// A discrete cash dividend: a known amount paid at a known future time.
// Added to a contract via OptionParams::addDividend, which validates it.
struct Dividend {
    double time;    // Payment time in years from now, in (0, T)
    double amount;  // Cash amount per share, non-negative
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
    std::vector<Dividend> discreteDividends;  // Cash dividends; empty by default
    std::vector<double> exerciseDates;        // Bermudan exercise dates; empty otherwise

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

    // Register a discrete cash dividend. Validated here so that engines can
    // rely on every stored dividend being payable strictly inside (0, T).
    void addDividend(double time, double amount) {
        if (time <= 0.0 || time >= timeToMaturity) {
            throw std::invalid_argument("Dividend time must lie strictly inside (0, maturity)");
        }
        if (amount < 0.0) {
            throw std::invalid_argument("Dividend amount must be non-negative");
        }
        discreteDividends.push_back({time, amount});
    }

    [[nodiscard]] bool hasDiscreteDividends() const noexcept {
        return !discreteDividends.empty();
    }

    // Register a Bermudan exercise date, validated into (0, T]. Exercise at
    // expiry itself is always available through the terminal payoff, so a
    // date equal to T is legal but redundant.
    void addExerciseDate(double time) {
        if (time <= 0.0 || time > timeToMaturity) {
            throw std::invalid_argument("Exercise date must lie inside (0, maturity]");
        }
        exerciseDates.push_back(time);
    }

    // Present value, discounted at the risk-free rate, of all dividends paid
    // after valuation time t (t = 0 gives the full escrowed adjustment).
    [[nodiscard]] double dividendPVAfter(double t) const noexcept {
        double pv = 0.0;
        for (const Dividend& div : discreteDividends) {
            if (div.time > t) {
                pv += div.amount * std::exp(-riskFreeRate * (div.time - t));
            }
        }
        return pv;
    }

    // The same contract with the spot reduced by the dividends' present value
    // and the dividend list cleared: the escrowed-dividend transformation
    // European engines price on directly.
    [[nodiscard]] OptionParams escrowed() const {
        OptionParams stripped = *this;
        stripped.spotPrice = spotPrice - dividendPVAfter(0.0);
        stripped.discreteDividends.clear();
        stripped.validate();  // A dividend PV >= spot leaves nothing to price
        return stripped;
    }

    // Helper methods
    [[nodiscard]] constexpr bool isCall() const noexcept { return optionType == OptionType::CALL; }
    [[nodiscard]] constexpr bool isPut() const noexcept { return optionType == OptionType::PUT; }
    [[nodiscard]] constexpr bool isEuropean() const noexcept { return exerciseType == ExerciseType::EUROPEAN; }
    [[nodiscard]] constexpr bool isAmerican() const noexcept { return exerciseType == ExerciseType::AMERICAN; }
    [[nodiscard]] constexpr bool isBermudan() const noexcept { return exerciseType == ExerciseType::BERMUDAN; }

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
