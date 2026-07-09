#include "options/BinomialTree.h"
#include "options/BlackScholes.h"
#include "options/ConvergenceAnalyzer.h"
#include "options/MonteCarlo.h"
#include "options/PerformanceBenchmark.h"

#include <cmath>
#include <iostream>

int main() {
    using namespace Options;

    const OptionParams call(100.0, 100.0, 0.05, 0.20, 1.0, 0.0,
                            OptionType::CALL, ExerciseType::EUROPEAN);
    const OptionParams put(100.0, 100.0, 0.05, 0.20, 1.0, 0.0,
                           OptionType::PUT, ExerciseType::EUROPEAN);

    BlackScholes bs;
    BinomialTree crr(2000);
    MonteCarlo mc(100000, DiscretizationScheme::EULER, VarianceReduction::BOTH);

    PricingEngine* const engines[] = {&bs, &crr, &mc};

    const double bsCall = bs.price(call).price;
    const double bsPut = bs.price(put).price;

    for (PricingEngine* engine : engines) {
        const PricingResult c = engine->price(call);
        const PricingResult p = engine->price(put);

        std::cout << engine->getName() << "\n"
                  << "  Call: " << c.price;
        if (c.standardError > 0.0) {
            std::cout << "  (SE " << c.standardError
                      << ", |err| " << std::fabs(c.price - bsCall) << " = "
                      << std::fabs(c.price - bsCall) / c.standardError << " SE)";
        }
        std::cout << "\n  Put:  " << p.price;
        if (p.standardError > 0.0) {
            std::cout << "  (SE " << p.standardError
                      << ", |err| " << std::fabs(p.price - bsPut) << " = "
                      << std::fabs(p.price - bsPut) / p.standardError << " SE)";
        }
        std::cout << "\n  Put-call parity (C - P): " << (c.price - p.price)
                  << "\n\n";
    }

    ConvergenceAnalyzer analyzer;
    const std::vector<size_t> binomialSteps = {10, 25, 50, 100, 250, 500, 1000, 2000};
    const std::vector<size_t> monteCarloPaths = {1000, 5000, 10000, 50000, 100000};

    const auto analyses =
        analyzer.compareAllMethods(call, binomialSteps, monteCarloPaths);
    for (const auto& analysis : analyses) {
        std::cout << analysis.toString() << "\n";
    }

    analyzer.exportToCSV(analyses, "results/convergence.csv");
    std::cout << "Convergence data written to results/convergence.csv\n\n";

    PerformanceBenchmark bench;
    bench.setNumRuns(10);
    const auto benchResults =
        bench.benchmarkAllMethods(call, {100, 1000}, {10000, 50000});
    std::cout << bench.createComparisonTable(benchResults);

    bench.exportToCSV(benchResults, "results/benchmark.csv");
    std::cout << "\nBenchmark data written to results/benchmark.csv\n";

    return 0;
}
