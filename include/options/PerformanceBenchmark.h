#ifndef PERFORMANCE_BENCHMARK_H
#define PERFORMANCE_BENCHMARK_H

#include "Option.h"
#include <vector>
#include <string>
#include <chrono>
#include <functional>

namespace Options {

// Structure to hold benchmark results
struct BenchmarkResult {
    std::string methodName;
    size_t iterations;
    double price;
    double avgTime;             // Average time per run (ms)
    double minTime;             // Minimum time (ms)
    double maxTime;             // Maximum time (ms)
    double stdDevTime;          // Standard deviation of time
    size_t memoryUsed;          // Memory used (bytes)
    double throughput;          // Iterations per second

    std::string toString() const;
};

// Performance Benchmark class
class PerformanceBenchmark {
public:
    PerformanceBenchmark() : numRuns_(10) {}

    // Set number of benchmark runs
    void setNumRuns(size_t runs) { numRuns_ = runs; }

    // Benchmark a specific method
    BenchmarkResult benchmark(
        const std::string& methodName,
        std::function<PricingResult()> pricingFunction
    );

    // Benchmark all methods and compare
    std::vector<BenchmarkResult> benchmarkAllMethods(
        const OptionParams& params,
        const std::vector<size_t>& binomialSteps,
        const std::vector<size_t>& monteCarloPaths
    );

    // Export benchmark results to CSV
    void exportToCSV(const std::vector<BenchmarkResult>& results,
                     const std::string& filename);

    // Create performance comparison table
    std::string createComparisonTable(const std::vector<BenchmarkResult>& results);

private:
    size_t numRuns_;

    // Calculate statistics
    static double calculateMean(const std::vector<double>& values);
    static double calculateStdDev(const std::vector<double>& values, double mean);
    static double calculateMin(const std::vector<double>& values);
    static double calculateMax(const std::vector<double>& values);
};

// Simple timer class
class Timer {
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}

    void reset() {
        start_ = std::chrono::high_resolution_clock::now();
    }

    double elapsed() const {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_);
        return duration.count() / 1000.0; // Return milliseconds
    }

private:
    std::chrono::high_resolution_clock::time_point start_;
};

} // namespace Options

#endif // PERFORMANCE_BENCHMARK_H
