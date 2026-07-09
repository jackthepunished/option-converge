#include "options/PerformanceBenchmark.h"
#include "options/BinomialTree.h"
#include "options/BlackScholes.h"
#include "options/MonteCarlo.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace Options {

std::string BenchmarkResult::toString() const {
    std::ostringstream out;
    out << methodName << "\n"
        << "  Price: " << price << "\n"
        << "  Runs: " << iterations << "\n"
        << "  Time (ms): avg " << avgTime << ", min " << minTime
        << ", max " << maxTime << ", stddev " << stdDevTime << "\n"
        << "  Throughput: " << throughput << " pricings/s\n"
        << "  Memory: " << memoryUsed << " bytes\n";
    return out.str();
}

double PerformanceBenchmark::calculateMean(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double v : values) {
        sum += v;
    }
    return sum / static_cast<double>(values.size());
}

double PerformanceBenchmark::calculateStdDev(const std::vector<double>& values,
                                             double mean) {
    if (values.size() < 2) {
        return 0.0;
    }
    double sumSq = 0.0;
    for (double v : values) {
        const double dev = v - mean;
        sumSq += dev * dev;
    }
    return std::sqrt(sumSq / static_cast<double>(values.size() - 1));
}

double PerformanceBenchmark::calculateMin(const std::vector<double>& values) {
    return values.empty() ? 0.0 : *std::min_element(values.begin(), values.end());
}

double PerformanceBenchmark::calculateMax(const std::vector<double>& values) {
    return values.empty() ? 0.0 : *std::max_element(values.begin(), values.end());
}

BenchmarkResult PerformanceBenchmark::benchmark(
    const std::string& methodName, std::function<PricingResult()> pricingFunction) {
    if (numRuns_ == 0) {
        throw std::invalid_argument("Benchmark requires at least one run");
    }

    // Warm-up run: first call pays for cache misses and lazy allocations and
    // would otherwise skew minTime/maxTime.
    PricingResult last = pricingFunction();

    // The timer resolves microseconds, so a single sub-microsecond call (e.g.
    // analytic Black-Scholes) would measure as zero. Batch calls until one
    // batch is comfortably measurable, then report per-call time.
    size_t batchSize = 1;
    constexpr size_t kMaxBatch = 1u << 20;
    constexpr double kMinBatchMs = 1.0;
    while (batchSize < kMaxBatch) {
        Timer probe;
        for (size_t i = 0; i < batchSize; ++i) {
            last = pricingFunction();
        }
        if (probe.elapsed() >= kMinBatchMs) {
            break;
        }
        batchSize *= 2;
    }

    std::vector<double> times;
    times.reserve(numRuns_);
    for (size_t i = 0; i < numRuns_; ++i) {
        Timer timer;
        for (size_t j = 0; j < batchSize; ++j) {
            last = pricingFunction();
        }
        times.push_back(timer.elapsed() / static_cast<double>(batchSize));
    }

    BenchmarkResult result;
    result.methodName = methodName;
    result.iterations = numRuns_;
    result.price = last.price;
    result.avgTime = calculateMean(times);
    result.minTime = calculateMin(times);
    result.maxTime = calculateMax(times);
    result.stdDevTime = calculateStdDev(times, result.avgTime);
    result.memoryUsed = last.memoryUsed;
    result.throughput = (result.avgTime > 0.0) ? 1000.0 / result.avgTime : 0.0;
    return result;
}

std::vector<BenchmarkResult> PerformanceBenchmark::benchmarkAllMethods(
    const OptionParams& params, const std::vector<size_t>& binomialSteps,
    const std::vector<size_t>& monteCarloPaths) {
    std::vector<BenchmarkResult> results;

    if (params.isEuropean()) {
        BlackScholes bs;
        results.push_back(
            benchmark("Black-Scholes", [&]() { return bs.price(params); }));
    }

    for (size_t steps : binomialSteps) {
        BinomialTree tree(steps);
        results.push_back(benchmark(tree.getName(), [&]() { return tree.price(params); }));
    }

    if (params.isEuropean()) {
        for (size_t paths : monteCarloPaths) {
            MonteCarlo mc(paths);
            results.push_back(benchmark(mc.getName(), [&]() {
                // Same seed each run so every repetition does identical work.
                mc.setSeed(42);
                return mc.price(params);
            }));
        }
    }

    return results;
}

std::string PerformanceBenchmark::createComparisonTable(
    const std::vector<BenchmarkResult>& results) {
    std::ostringstream out;
    out << std::left << std::setw(52) << "Method" << std::right
        << std::setw(12) << "Price" << std::setw(12) << "Avg (ms)"
        << std::setw(12) << "Min (ms)" << std::setw(12) << "Max (ms)"
        << std::setw(16) << "Pricings/s" << "\n";
    out << std::string(116, '-') << "\n";

    out << std::fixed << std::setprecision(6);
    for (const BenchmarkResult& r : results) {
        out << std::left << std::setw(52) << r.methodName << std::right
            << std::setw(12) << r.price << std::setw(12) << r.avgTime
            << std::setw(12) << r.minTime << std::setw(12) << r.maxTime
            << std::setw(16) << std::setprecision(1) << r.throughput
            << std::setprecision(6) << "\n";
    }
    return out.str();
}

void PerformanceBenchmark::exportToCSV(const std::vector<BenchmarkResult>& results,
                                       const std::string& filename) {
    const std::filesystem::path path(filename);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }

    out << "method,runs,price,avg_time_ms,min_time_ms,max_time_ms,"
           "stddev_time_ms,throughput_per_s,memory_bytes\n";
    for (const BenchmarkResult& r : results) {
        out << '"' << r.methodName << '"' << ',' << r.iterations << ','
            << r.price << ',' << r.avgTime << ',' << r.minTime << ','
            << r.maxTime << ',' << r.stdDevTime << ',' << r.throughput << ','
            << r.memoryUsed << '\n';
    }
}

} // namespace Options
