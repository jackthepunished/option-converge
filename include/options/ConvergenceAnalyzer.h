#ifndef CONVERGENCE_ANALYZER_H
#define CONVERGENCE_ANALYZER_H

#include "Option.h"
#include "PricingEngine.h"
#include <vector>
#include <string>
#include <memory>

namespace Options {

// Structure to hold convergence data point
struct ConvergencePoint {
    size_t iterations;          // Number of steps/paths
    double price;               // Calculated price
    double error;               // Error vs reference price
    double relativeError;       // Relative error (%)
    double computationTime;     // Time in milliseconds
    size_t memoryUsed;          // Memory in bytes
    Greeks greeks;              // Greeks at this point

    ConvergencePoint() : iterations(0), price(0), error(0),
                        relativeError(0), computationTime(0), memoryUsed(0) {}
};

// Structure to hold full convergence analysis
struct ConvergenceAnalysis {
    std::string methodName;
    std::vector<ConvergencePoint> points;
    double referencePrice;      // Black-Scholes reference
    double rmse;                // Root Mean Square Error
    double avgComputationTime;  // Average computation time
    size_t convergenceIteration; // Iteration where error < 0.01%

    std::string toString() const;
};

// Convergence Analyzer class
class ConvergenceAnalyzer {
public:
    ConvergenceAnalyzer() = default;

    // Analyze convergence for Binomial Tree
    ConvergenceAnalysis analyzeBinomialConvergence(
        const OptionParams& params,
        const std::vector<size_t>& stepSizes,
        double referencePrice
    );

    // Analyze convergence for Monte Carlo
    ConvergenceAnalysis analyzeMonteCarloConvergence(
        const OptionParams& params,
        const std::vector<size_t>& pathCounts,
        double referencePrice,
        DiscretizationScheme scheme = DiscretizationScheme::EULER,
        VarianceReduction varRed = VarianceReduction::ANTITHETIC
    );

    // Compare all methods
    std::vector<ConvergenceAnalysis> compareAllMethods(
        const OptionParams& params,
        const std::vector<size_t>& binomialSteps,
        const std::vector<size_t>& monteCarloPaths
    );

    // Export results to CSV
    void exportToCSV(const std::vector<ConvergenceAnalysis>& analyses,
                     const std::string& filename);

    // Calculate RMSE
    static double calculateRMSE(const std::vector<ConvergencePoint>& points,
                               double referencePrice);

    // Find convergence iteration (error < threshold)
    static size_t findConvergenceIteration(const std::vector<ConvergencePoint>& points,
                                          double threshold = 0.0001);

private:
    // Helper to calculate errors
    void calculateErrors(std::vector<ConvergencePoint>& points,
                        double referencePrice);
};

} // namespace Options

#endif // CONVERGENCE_ANALYZER_H
