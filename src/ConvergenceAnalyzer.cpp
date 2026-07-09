#include "options/ConvergenceAnalyzer.h"
#include "options/BinomialTree.h"
#include "options/BlackScholes.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace Options {

std::string ConvergenceAnalysis::toString() const {
    std::ostringstream out;
    out << "Convergence Analysis: " << methodName << "\n"
        << "  Reference price: " << referencePrice << "\n"
        << "  RMSE: " << rmse << "\n"
        << "  Avg computation time: " << avgComputationTime << " ms\n";
    if (convergenceIteration > 0) {
        out << "  Converged (|rel err| < 0.01%) at: " << convergenceIteration
            << " iterations\n";
    } else {
        out << "  Did not reach 0.01% relative error\n";
    }
    out << "  Points:\n";
    for (const ConvergencePoint& p : points) {
        out << "    n=" << p.iterations << "  price=" << p.price
            << "  err=" << p.error << "  relErr=" << p.relativeError << "%"
            << "  time=" << p.computationTime << "ms\n";
    }
    return out.str();
}

void ConvergenceAnalyzer::calculateErrors(std::vector<ConvergencePoint>& points,
                                          double referencePrice) {
    for (ConvergencePoint& p : points) {
        p.error = p.price - referencePrice;
        p.relativeError = (referencePrice != 0.0)
                              ? 100.0 * std::fabs(p.error) / std::fabs(referencePrice)
                              : 0.0;
    }
}

double ConvergenceAnalyzer::calculateRMSE(const std::vector<ConvergencePoint>& points,
                                          double referencePrice) {
    if (points.empty()) {
        return 0.0;
    }
    double sumSq = 0.0;
    for (const ConvergencePoint& p : points) {
        const double err = p.price - referencePrice;
        sumSq += err * err;
    }
    return std::sqrt(sumSq / static_cast<double>(points.size()));
}

size_t ConvergenceAnalyzer::findConvergenceIteration(
    const std::vector<ConvergencePoint>& points, double threshold) {
    // First point from which every later point also stays inside the
    // threshold — a single lucky crossing does not count as convergence.
    for (size_t i = 0; i < points.size(); ++i) {
        bool stable = true;
        for (size_t j = i; j < points.size(); ++j) {
            if (points[j].relativeError > 100.0 * threshold) {
                stable = false;
                break;
            }
        }
        if (stable) {
            return points[i].iterations;
        }
    }
    return 0;
}

ConvergenceAnalysis ConvergenceAnalyzer::analyzeBinomialConvergence(
    const OptionParams& params, const std::vector<size_t>& stepSizes,
    double referencePrice) {
    ConvergenceAnalysis analysis;
    analysis.referencePrice = referencePrice;

    BinomialTree engine(stepSizes.empty() ? 100 : stepSizes.front());
    analysis.methodName = "Binomial Tree (CRR)";

    double totalTime = 0.0;
    for (size_t steps : stepSizes) {
        engine.setSteps(steps);
        const PricingResult result = engine.price(params);

        ConvergencePoint point;
        point.iterations = steps;
        point.price = result.price;
        point.computationTime = result.computationTime;
        point.memoryUsed = result.memoryUsed;
        point.greeks = result.greeks;

        totalTime += result.computationTime;
        analysis.points.push_back(point);
    }

    calculateErrors(analysis.points, referencePrice);
    analysis.rmse = calculateRMSE(analysis.points, referencePrice);
    analysis.avgComputationTime =
        analysis.points.empty() ? 0.0
                                : totalTime / static_cast<double>(analysis.points.size());
    analysis.convergenceIteration = findConvergenceIteration(analysis.points);

    return analysis;
}

ConvergenceAnalysis ConvergenceAnalyzer::analyzeMonteCarloConvergence(
    const OptionParams& params, const std::vector<size_t>& pathCounts,
    double referencePrice, DiscretizationScheme scheme, VarianceReduction varRed) {
    ConvergenceAnalysis analysis;
    analysis.referencePrice = referencePrice;

    MonteCarlo engine(pathCounts.empty() ? 10000 : pathCounts.front(), scheme, varRed);

    // Engine name without the path count: the sweep varies it.
    std::string label = "Monte Carlo (";
    label += (scheme == DiscretizationScheme::EULER) ? "Euler" : "Milstein";
    switch (varRed) {
        case VarianceReduction::NONE: label += ")"; break;
        case VarianceReduction::ANTITHETIC: label += ", antithetic)"; break;
        case VarianceReduction::CONTROL_VARIATE: label += ", control variate)"; break;
        case VarianceReduction::BOTH: label += ", antithetic + control variate)"; break;
    }
    analysis.methodName = label;

    double totalTime = 0.0;
    for (size_t paths : pathCounts) {
        engine.setNumPaths(paths);
        // Same seed for every sample size: the error series then reflects the
        // effect of adding paths, not a fresh random draw each time.
        engine.setSeed(42);
        const PricingResult result = engine.price(params);

        ConvergencePoint point;
        point.iterations = paths;
        point.price = result.price;
        point.computationTime = result.computationTime;
        point.memoryUsed = result.memoryUsed;
        point.greeks = result.greeks;

        totalTime += result.computationTime;
        analysis.points.push_back(point);
    }

    calculateErrors(analysis.points, referencePrice);
    analysis.rmse = calculateRMSE(analysis.points, referencePrice);
    analysis.avgComputationTime =
        analysis.points.empty() ? 0.0
                                : totalTime / static_cast<double>(analysis.points.size());
    analysis.convergenceIteration = findConvergenceIteration(analysis.points);

    return analysis;
}

std::vector<ConvergenceAnalysis> ConvergenceAnalyzer::compareAllMethods(
    const OptionParams& params, const std::vector<size_t>& binomialSteps,
    const std::vector<size_t>& monteCarloPaths) {
    if (!params.isEuropean()) {
        throw std::invalid_argument(
            "compareAllMethods needs a European option: the Black-Scholes "
            "reference price is only defined for European exercise.");
    }

    BlackScholes reference;
    const double referencePrice = reference.price(params).price;

    std::vector<ConvergenceAnalysis> analyses;
    analyses.push_back(
        analyzeBinomialConvergence(params, binomialSteps, referencePrice));
    analyses.push_back(analyzeMonteCarloConvergence(
        params, monteCarloPaths, referencePrice));
    return analyses;
}

void ConvergenceAnalyzer::exportToCSV(
    const std::vector<ConvergenceAnalysis>& analyses, const std::string& filename) {
    const std::filesystem::path path(filename);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }

    out << "method,iterations,price,reference_price,error,relative_error_pct,"
           "computation_time_ms,memory_bytes\n";
    for (const ConvergenceAnalysis& analysis : analyses) {
        for (const ConvergencePoint& p : analysis.points) {
            out << '"' << analysis.methodName << '"' << ','
                << p.iterations << ',' << p.price << ','
                << analysis.referencePrice << ',' << p.error << ','
                << p.relativeError << ',' << p.computationTime << ','
                << p.memoryUsed << '\n';
        }
    }
}

} // namespace Options
