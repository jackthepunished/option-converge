// Test suite for the options pricing library. Hand-rolled assertions keep the
// project dependency-free; the process exit code is the number of failures.

#include "options/BinomialTree.h"
#include "options/BlackScholes.h"
#include "options/ConvergenceAnalyzer.h"
#include "options/MonteCarlo.h"
#include "options/PerformanceBenchmark.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

int failures = 0;
int checks = 0;

void check(bool condition, const std::string& label) {
    ++checks;
    if (!condition) {
        ++failures;
        std::cout << "FAIL: " << label << "\n";
    }
}

void checkNear(double actual, double expected, double tolerance,
               const std::string& label) {
    ++checks;
    if (std::fabs(actual - expected) > tolerance) {
        ++failures;
        std::cout << "FAIL: " << label << "  expected " << expected << " +/- "
                  << tolerance << ", got " << actual << "\n";
    }
}

template <typename Fn>
void checkThrows(Fn fn, const std::string& label) {
    ++checks;
    bool threw = false;
    try {
        fn();
    } catch (const std::exception&) {
        threw = true;
    }
    if (!threw) {
        ++failures;
        std::cout << "FAIL: " << label << "  (no exception thrown)\n";
    }
}

// Reference option: S=100, K=100, r=5%, sigma=20%, T=1, no dividends.
const Options::OptionParams kCall(100.0, 100.0, 0.05, 0.20, 1.0, 0.0,
                                  Options::OptionType::CALL,
                                  Options::ExerciseType::EUROPEAN);
const Options::OptionParams kPut(100.0, 100.0, 0.05, 0.20, 1.0, 0.0,
                                 Options::OptionType::PUT,
                                 Options::ExerciseType::EUROPEAN);

constexpr double kBsCall = 10.450584;
constexpr double kBsPut = 5.573526;
const double kParity = 100.0 - 100.0 * std::exp(-0.05);  // S - K*exp(-rT)

void testOptionParamsValidation() {
    using namespace Options;
    checkThrows([] { OptionParams(-1.0, 100.0); }, "negative spot rejected");
    checkThrows([] { OptionParams(100.0, 0.0); }, "zero strike rejected");
    checkThrows([] { OptionParams(100.0, 100.0, 0.05, -0.1); },
                "negative volatility rejected");
    checkThrows([] { OptionParams(100.0, 100.0, 0.05, 0.2, 0.0); },
                "zero maturity rejected");
}

void testBlackScholes() {
    using namespace Options;
    BlackScholes bs;

    const double call = bs.price(kCall).price;
    const double put = bs.price(kPut).price;
    checkNear(call, kBsCall, 1e-5, "BS call price");
    checkNear(put, kBsPut, 1e-5, "BS put price");
    checkNear(call - put, kParity, 1e-9, "BS put-call parity");

    // Textbook Greeks for the reference option.
    const Greeks callGreeks = bs.calculateGreeks(kCall);
    checkNear(callGreeks.delta, 0.636831, 1e-3, "BS call delta");
    checkNear(callGreeks.gamma, 0.018762, 1e-3, "BS gamma");
    checkNear(callGreeks.vega, 0.375240, 1e-2, "BS vega");
    checkNear(callGreeks.theta, -0.017573, 1e-4, "BS call theta (daily)");
    checkNear(callGreeks.rho, 0.532325, 1e-3, "BS call rho (per 1%)");

    const Greeks putGreeks = bs.calculateGreeks(kPut);
    checkNear(putGreeks.delta, -0.363169, 1e-3, "BS put delta");
    checkNear(putGreeks.theta, -0.004542, 1e-4, "BS put theta (daily)");

    // Finite-difference Greeks from the lattice must agree with the analytic
    // values under the same conventions (vega/rho per 1%, theta per day).
    BinomialTree tree(1000);
    const Greeks fdGreeks = tree.calculateGreeks(kCall);
    checkNear(fdGreeks.delta, callGreeks.delta, 1e-3, "CRR delta matches BS");
    checkNear(fdGreeks.vega, callGreeks.vega, 1e-2, "CRR vega matches BS");
    checkNear(fdGreeks.theta, callGreeks.theta, 1e-3, "CRR theta matches BS");
    checkNear(fdGreeks.rho, callGreeks.rho, 1e-2, "CRR rho matches BS");

    OptionParams american = kCall;
    american.exerciseType = ExerciseType::AMERICAN;
    checkThrows([&] { (void)BlackScholes().price(american); },
                "BS rejects American options");
}

void testBinomialTree() {
    using namespace Options;
    BinomialTree tree(1000);

    const double call = tree.price(kCall).price;
    const double put = tree.price(kPut).price;
    checkNear(call, kBsCall, 0.01, "CRR(1000) call converges to BS");
    checkNear(put, kBsPut, 0.01, "CRR(1000) put converges to BS");
    checkNear(call - put, kParity, 1e-6, "CRR put-call parity");

    // Without dividends, early exercise of an American call is never optimal,
    // so its lattice price must match the European call. An American put
    // carries a strictly positive early-exercise premium.
    OptionParams americanCall = kCall;
    americanCall.exerciseType = ExerciseType::AMERICAN;
    checkNear(tree.price(americanCall).price, call, 1e-9,
              "American call equals European call (q=0)");

    OptionParams americanPut = kPut;
    americanPut.exerciseType = ExerciseType::AMERICAN;
    const double amPut = tree.price(americanPut).price;
    check(amPut > put + 0.05, "American put has early-exercise premium");
    check(amPut < put + 2.0, "American put premium is plausible");

    // Deep in-the-money American put: immediate exercise dominates, so the
    // price must be at least intrinsic value.
    OptionParams deepPut(50.0, 100.0, 0.05, 0.20, 1.0, 0.0, OptionType::PUT,
                         ExerciseType::AMERICAN);
    check(tree.price(deepPut).price >= 50.0 - 1e-9,
          "deep ITM American put >= intrinsic");

    checkThrows([] { BinomialTree(0); }, "zero steps rejected");
}

void testMonteCarlo() {
    using namespace Options;

    const VarianceReduction modes[] = {
        VarianceReduction::NONE, VarianceReduction::ANTITHETIC,
        VarianceReduction::CONTROL_VARIATE, VarianceReduction::BOTH};
    const char* labels[] = {"none", "antithetic", "control variate", "both"};

    for (int i = 0; i < 4; ++i) {
        MonteCarlo mc(50000, DiscretizationScheme::EULER, modes[i]);
        const PricingResult call = mc.price(kCall);
        check(call.standardError > 0.0,
              std::string("MC SE positive (") + labels[i] + ")");
        checkNear(call.price, kBsCall, 3.0 * call.standardError,
                  std::string("MC call within 3 SE (") + labels[i] + ")");
    }

    // Milstein scheme must agree too.
    MonteCarlo milstein(50000, DiscretizationScheme::MILSTEIN,
                        VarianceReduction::ANTITHETIC);
    const PricingResult mCall = milstein.price(kCall);
    checkNear(mCall.price, kBsCall, 3.0 * mCall.standardError,
              "MC Milstein call within 3 SE");

    // Variance reduction has to actually reduce variance.
    MonteCarlo basic(50000, DiscretizationScheme::EULER, VarianceReduction::NONE);
    MonteCarlo both(50000, DiscretizationScheme::EULER, VarianceReduction::BOTH);
    const double seBasic = basic.price(kCall).standardError;
    const double seBoth = both.price(kCall).standardError;
    check(seBoth < seBasic,
          "antithetic + control variate reduces standard error");

    // Same seed, same estimate: pricing must be reproducible.
    MonteCarlo a(20000), b(20000);
    checkNear(a.price(kCall).price, b.price(kCall).price, 1e-12,
              "same seed gives identical MC price");

    MonteCarlo mc(10000);
    OptionParams american = kPut;
    american.exerciseType = ExerciseType::AMERICAN;
    checkThrows([&] { (void)mc.price(american); }, "MC rejects American options");
}

void testConvergenceAnalyzer() {
    using namespace Options;
    ConvergenceAnalyzer analyzer;

    const auto binomial = analyzer.analyzeBinomialConvergence(
        kCall, {10, 100, 1000}, kBsCall);
    check(binomial.points.size() == 3, "binomial analysis has all points");
    check(std::fabs(binomial.points.back().error) <
              std::fabs(binomial.points.front().error),
          "binomial error shrinks with more steps");
    check(binomial.rmse > 0.0, "binomial RMSE positive");

    const auto mc = analyzer.analyzeMonteCarloConvergence(
        kCall, {1000, 10000}, kBsCall);
    check(mc.points.size() == 2, "MC analysis has all points");

    // RMSE and convergence-iteration helpers on hand-built data.
    std::vector<ConvergencePoint> points(2);
    points[0].iterations = 10;
    points[0].price = kBsCall + 0.3;
    points[1].iterations = 20;
    points[1].price = kBsCall - 0.4;
    checkNear(ConvergenceAnalyzer::calculateRMSE(points, kBsCall),
              std::sqrt((0.09 + 0.16) / 2.0), 1e-12, "RMSE formula");

    points[0].relativeError = 0.5;   // percent
    points[1].relativeError = 0.001; // percent, below 0.01% threshold
    check(ConvergenceAnalyzer::findConvergenceIteration(points) == 20,
          "convergence iteration found");
    points[1].relativeError = 0.5;
    check(ConvergenceAnalyzer::findConvergenceIteration(points) == 0,
          "no convergence reported when never under threshold");

    checkThrows(
        [&] {
            OptionParams american = kCall;
            american.exerciseType = ExerciseType::AMERICAN;
            (void)analyzer.compareAllMethods(american, {10}, {1000});
        },
        "compareAllMethods rejects American reference");
}

void testPerformanceBenchmark() {
    using namespace Options;
    PerformanceBenchmark bench;
    bench.setNumRuns(3);

    BlackScholes bs;
    const BenchmarkResult result =
        bench.benchmark("Black-Scholes", [&] { return bs.price(kCall); });
    checkNear(result.price, kBsCall, 1e-5, "benchmark reports correct price");
    check(result.avgTime > 0.0, "benchmark measures nonzero time");
    check(result.throughput > 0.0, "benchmark reports throughput");
    check(result.minTime <= result.avgTime && result.avgTime <= result.maxTime,
          "benchmark stats ordered");
}

} // namespace

int main() {
    testOptionParamsValidation();
    testBlackScholes();
    testBinomialTree();
    testMonteCarlo();
    testConvergenceAnalyzer();
    testPerformanceBenchmark();

    std::cout << (checks - failures) << "/" << checks << " checks passed\n";
    return failures;
}
