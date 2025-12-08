#include "options/Option.h"
#include <sstream>
#include <iomanip>

namespace Options {

std::string OptionParams::toString() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4);

    oss << "Option Parameters:\n";
    oss << "  Type: " << (isCall() ? "Call" : "Put") << "\n";
    oss << "  Exercise: " << (isEuropean() ? "European" : "American") << "\n";
    oss << "  Spot Price (S): " << spotPrice << "\n";
    oss << "  Strike Price (K): " << strikePrice << "\n";
    oss << "  Risk-Free Rate (r): " << (riskFreeRate * 100) << "%\n";
    oss << "  Volatility (σ): " << (volatility * 100) << "%\n";
    oss << "  Time to Maturity (T): " << timeToMaturity << " years\n";
    oss << "  Dividend Yield (q): " << (dividendYield * 100) << "%";

    return oss.str();
}

std::string Greeks::toString() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);

    oss << "Greeks:\n";
    oss << "  Delta (Δ): " << delta << "\n";
    oss << "  Gamma (Γ): " << gamma << "\n";
    oss << "  Vega (ν): " << vega << "\n";
    oss << "  Theta (Θ): " << theta << "\n";
    oss << "  Rho (ρ): " << rho;

    return oss.str();
}

std::string PricingResult::toString() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);

    oss << "Pricing Result:\n";
    oss << "  Price: " << price << "\n";

    if (standardError > 0) {
        oss << "  Standard Error: " << standardError << "\n";
    }

    oss << "  Computation Time: " << std::setprecision(2)
        << computationTime << " ms\n";

    if (memoryUsed > 0) {
        double memoryMB = memoryUsed / (1024.0 * 1024.0);
        oss << "  Memory Used: " << std::setprecision(2)
            << memoryMB << " MB\n";
    }

    oss << greeks.toString();

    return oss.str();
}

} // namespace Options
