// Test suite for the options pricing library. Hand-rolled assertions keep the
// project dependency-free; the process exit code is the number of failures.

#include "options/AsianOption.h"
#include "options/BarrierOption.h"
#include "options/BinomialTree.h"
#include "options/BlackScholes.h"
#include "options/ConvergenceAnalyzer.h"
#include "options/FiniteDifference.h"
#include "options/ImpliedVolatility.h"
#include "options/LongstaffSchwartz.h"
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

void testLongstaffSchwartz() {
    using namespace Options;
    BinomialTree tree(2000);
    LongstaffSchwartz lsmc(50000, 50);

    OptionParams americanPut = kPut;
    americanPut.exerciseType = ExerciseType::AMERICAN;
    OptionParams americanCall = kCall;
    americanCall.exerciseType = ExerciseType::AMERICAN;

    // The regression estimate must land on the lattice's American put value:
    // two unrelated treatments of early exercise, one price. LSMC carries a
    // small low bias from the finite basis, hence the additive allowance.
    const auto putEst = lsmc.price(americanPut);
    const double treePut = tree.price(americanPut).price;
    check(putEst.standardError > 0.0, "LSMC reports standard error");
    checkNear(putEst.price, treePut, 3.0 * putEst.standardError + 0.05,
              "LSMC American put matches lattice");

    // The early-exercise premium over the European put must survive.
    check(putEst.price > kBsPut + 0.3,
          "LSMC American put carries early-exercise premium");

    // Without dividends, early exercise of a call is never optimal, so the
    // LSMC American call must reproduce the European Black-Scholes value.
    const auto callEst = lsmc.price(americanCall);
    checkNear(callEst.price, kBsCall, 3.0 * callEst.standardError + 0.05,
              "LSMC American call equals European (q=0)");

    // With a heavy dividend yield early exercise of the call has value, and
    // LSMC must agree with the lattice about how much.
    OptionParams divCall(100.0, 100.0, 0.05, 0.20, 1.0, 0.08, OptionType::CALL,
                         ExerciseType::AMERICAN);
    const auto divEst = lsmc.price(divCall);
    checkNear(divEst.price, tree.price(divCall).price,
              3.0 * divEst.standardError + 0.05,
              "LSMC dividend-paying American call matches lattice");

    // Deep in the money, immediate exercise dominates and the intrinsic
    // floor at t=0 must hold exactly.
    OptionParams deepPut(50.0, 100.0, 0.05, 0.20, 1.0, 0.0, OptionType::PUT,
                         ExerciseType::AMERICAN);
    check(lsmc.price(deepPut).price >= 50.0 - 1e-12,
          "LSMC deep ITM American put >= intrinsic");

    // Same seed, same estimate.
    LongstaffSchwartz a(20000, 25), b(20000, 25);
    checkNear(a.price(americanPut).price, b.price(americanPut).price, 1e-12,
              "same seed gives identical LSMC price");

    checkThrows([&] { (void)lsmc.price(kPut); }, "LSMC rejects European options");
    checkThrows([] { LongstaffSchwartz(1, 50); }, "degenerate path count rejected");
    checkThrows([] { LongstaffSchwartz(1000, 0); }, "zero exercise dates rejected");
}

void testAsianOptions() {
    using namespace Options;
    BlackScholes bs;

    // With a single fixing the average IS the terminal spot, so the discrete
    // geometric closed form must reproduce Black-Scholes exactly.
    AsianParams oneFixCall(kCall, AveragingType::GEOMETRIC, 1);
    AsianParams oneFixPut(kPut, AveragingType::GEOMETRIC, 1);
    checkNear(AnalyticGeometricAsian::price(oneFixCall), kBsCall, 1e-6,
              "geometric Asian with one fixing is the vanilla call");
    checkNear(AnalyticGeometricAsian::price(oneFixPut), kBsPut, 1e-6,
              "geometric Asian with one fixing is the vanilla put");

    // Averaging damps volatility, so the Asian call is worth less than the
    // vanilla, and its value falls as more fixings enter the average.
    AsianParams geo12(kCall, AveragingType::GEOMETRIC, 12);
    AsianParams geo50(kCall, AveragingType::GEOMETRIC, 50);
    const double geo12Price = AnalyticGeometricAsian::price(geo12);
    const double geo50Price = AnalyticGeometricAsian::price(geo50);
    check(geo12Price < kBsCall, "geometric Asian call below vanilla call");
    check(geo50Price < geo12Price, "more fixings, lower Asian call value");

    // Put-call parity in the geometric measure: C - P = disc * (E[G] - K).
    AsianParams geo12Put(kPut, AveragingType::GEOMETRIC, 12);
    const double geoPutPrice = AnalyticGeometricAsian::price(geo12Put);
    {
        // Recompute E[G] from the same moments the formula uses.
        const double nd = 12.0;
        const double m = std::log(100.0) + (0.05 - 0.5 * 0.04) * 1.0 * (nd + 1.0) / (2.0 * nd);
        const double v2 = 0.04 * 1.0 * (nd + 1.0) * (2.0 * nd + 1.0) / (6.0 * nd * nd);
        const double parity = std::exp(-0.05) * (std::exp(m + 0.5 * v2) - 100.0);
        checkNear(geo12Price - geoPutPrice, parity, 1e-10, "geometric Asian put-call parity");
    }

    // Monte Carlo geometric pricing must agree with its own closed form.
    AsianMonteCarlo mc(100000);
    {
        const auto est = mc.price(geo12);
        check(est.standardError > 0.0, "Asian MC reports standard error");
        checkNear(est.price, geo12Price, 3.0 * est.standardError,
                  "MC geometric Asian within 3 SE of closed form");
    }

    // The arithmetic mean dominates the geometric (AM-GM), so the arithmetic
    // call must be worth at least the geometric call.
    AsianParams arith12(kCall, AveragingType::ARITHMETIC, 12);
    const auto arithEst = mc.price(arith12);
    check(arithEst.price > geo12Price,
          "arithmetic Asian call above geometric (AM-GM)");
    check(arithEst.price < kBsCall, "arithmetic Asian call below vanilla call");

    // The geometric control variate has to earn its keep.
    const auto plain = mc.price(arith12, false);
    check(arithEst.standardError < 0.2 * plain.standardError,
          "geometric control variate slashes arithmetic SE");
    checkNear(arithEst.price, plain.price, 3.0 * plain.standardError,
              "CV and plain estimates agree");

    // Same seed, same estimate.
    AsianMonteCarlo mcA(20000), mcB(20000);
    checkNear(mcA.price(arith12).price, mcB.price(arith12).price, 1e-12,
              "same seed gives identical Asian MC price");

    checkThrows([] { AsianParams(kCall, AveragingType::ARITHMETIC, 0); },
                "zero fixings rejected");
    checkThrows(
        [] {
            OptionParams american = kCall;
            american.exerciseType = ExerciseType::AMERICAN;
            AsianParams(american, AveragingType::GEOMETRIC, 12);
        },
        "American Asian rejected");
    checkThrows(
        [] {
            AsianParams arith(kCall, AveragingType::ARITHMETIC, 12);
            (void)AnalyticGeometricAsian::price(arith);
        },
        "closed form refuses arithmetic averaging");
}

void testBarrierOptions() {
    using namespace Options;
    BlackScholes bs;

    // Exact structural identities. An up-and-in call with the barrier at or
    // below the strike IS the vanilla call: any path ending in the money has
    // necessarily crossed the barrier. Its knock-out complement is worthless.
    {
        OptionParams call110 = kCall;
        call110.strikePrice = 110.0;
        const double vanilla = bs.price(call110).price;
        BarrierParams upIn(call110, BarrierType::UP_AND_IN, 110.0);
        BarrierParams upOut(call110, BarrierType::UP_AND_OUT, 110.0);
        checkNear(AnalyticBarrier::price(upIn), vanilla, 1e-10,
                  "up-in call with B<=K equals vanilla");
        checkNear(AnalyticBarrier::price(upOut), 0.0, 1e-10,
                  "up-out call with B<=K is worthless");
    }
    {
        OptionParams put85 = kPut;
        put85.strikePrice = 85.0;
        const double vanilla = bs.price(put85).price;
        BarrierParams downIn(put85, BarrierType::DOWN_AND_IN, 90.0);
        BarrierParams downOut(put85, BarrierType::DOWN_AND_OUT, 90.0);
        checkNear(AnalyticBarrier::price(downIn), vanilla, 1e-10,
                  "down-in put with B>=K equals vanilla");
        checkNear(AnalyticBarrier::price(downOut), 0.0, 1e-10,
                  "down-out put with B>=K is worthless");
    }

    // A barrier the process cannot plausibly reach changes nothing.
    BarrierParams farOut(kCall, BarrierType::DOWN_AND_OUT, 1.0);
    BarrierParams farIn(kCall, BarrierType::DOWN_AND_IN, 1.0);
    checkNear(AnalyticBarrier::price(farOut), kBsCall, 1e-6,
              "far barrier knock-out equals vanilla");
    checkNear(AnalyticBarrier::price(farIn), 0.0, 1e-6,
              "far barrier knock-in is worthless");

    // Tightening a knock-out barrier can only destroy value.
    const double out80 =
        AnalyticBarrier::price(BarrierParams(kCall, BarrierType::DOWN_AND_OUT, 80.0));
    const double out95 =
        AnalyticBarrier::price(BarrierParams(kCall, BarrierType::DOWN_AND_OUT, 95.0));
    check(out95 < out80 && out80 < kBsCall + 1e-9,
          "down-out call value decreases as barrier rises");

    // Monte Carlo with the Brownian-bridge correction must agree with the
    // continuously monitored closed form — two unrelated methods, one price.
    BarrierMonteCarlo mc(100000, 64);
    {
        BarrierParams downOut(kCall, BarrierType::DOWN_AND_OUT, 90.0);
        const auto est = mc.price(downOut);
        check(est.standardError > 0.0, "barrier MC reports standard error");
        checkNear(est.price, AnalyticBarrier::price(downOut),
                  3.0 * est.standardError + 0.01, "MC down-out call vs analytic");
    }
    {
        BarrierParams upOut(kCall, BarrierType::UP_AND_OUT, 120.0);
        const auto est = mc.price(upOut);
        checkNear(est.price, AnalyticBarrier::price(upOut),
                  3.0 * est.standardError + 0.01, "MC up-out call vs analytic");
    }
    {
        BarrierParams downIn(kPut, BarrierType::DOWN_AND_IN, 90.0);
        const auto est = mc.price(downIn);
        checkNear(est.price, AnalyticBarrier::price(downIn),
                  3.0 * est.standardError + 0.01, "MC down-in put vs analytic");
    }

    // In-out parity holds path by path under the same seed, so the sum is
    // exactly the plain MC vanilla estimate and must sit near Black-Scholes.
    const auto mcIn = mc.price(BarrierParams(kCall, BarrierType::DOWN_AND_IN, 90.0));
    const auto mcOut = mc.price(BarrierParams(kCall, BarrierType::DOWN_AND_OUT, 90.0));
    checkNear(mcIn.price + mcOut.price, kBsCall, 0.15, "MC in + out = vanilla");

    // Already-breached contracts degenerate immediately.
    BarrierParams breachedOut(kCall, BarrierType::DOWN_AND_OUT, 100.0);
    BarrierParams breachedIn(kCall, BarrierType::UP_AND_IN, 95.0);
    checkNear(AnalyticBarrier::price(breachedOut), 0.0, 1e-12,
              "breached knock-out is worthless");
    checkNear(AnalyticBarrier::price(breachedIn), kBsCall, 1e-5,
              "breached knock-in is the vanilla option");
    checkNear(mc.price(breachedOut).price, 0.0, 1e-12, "MC breached knock-out is zero");

    checkThrows([] { BarrierParams(kCall, BarrierType::UP_AND_OUT, -5.0); },
                "non-positive barrier rejected");
    checkThrows(
        [] {
            OptionParams american = kCall;
            american.exerciseType = ExerciseType::AMERICAN;
            BarrierParams(american, BarrierType::UP_AND_OUT, 120.0);
        },
        "American barrier rejected");
}

void testFiniteDifference() {
    using namespace Options;

    // All three schemes must converge to the analytic European price. The
    // explicit scheme refines its own time axis to stay inside its stability
    // bound, so the same grid request is safe for every scheme.
    FiniteDifference explicitFd(200, 200, FDScheme::EXPLICIT);
    FiniteDifference implicitFd(200, 200, FDScheme::IMPLICIT);
    FiniteDifference cn(200, 200, FDScheme::CRANK_NICOLSON);

    checkNear(explicitFd.price(kCall).price, kBsCall, 0.01, "FD explicit call vs BS");
    checkNear(implicitFd.price(kCall).price, kBsCall, 0.01, "FD implicit call vs BS");
    checkNear(cn.price(kCall).price, kBsCall, 0.005, "FD Crank-Nicolson call vs BS");
    checkNear(cn.price(kPut).price, kBsPut, 0.005, "FD Crank-Nicolson put vs BS");

    checkNear(cn.price(kCall).price - cn.price(kPut).price, kParity, 5e-3,
              "FD put-call parity");

    // American exercise must agree with the lattice, the established
    // early-exercise reference in this library.
    BinomialTree tree(1000);
    OptionParams americanPut = kPut;
    americanPut.exerciseType = ExerciseType::AMERICAN;
    const double treeAmPut = tree.price(americanPut).price;
    checkNear(cn.price(americanPut).price, treeAmPut, 0.02,
              "FD American put matches lattice");
    check(cn.price(americanPut).price > cn.price(kPut).price + 0.05,
          "FD American put has early-exercise premium");

    // Without dividends early exercise of a call is never optimal.
    OptionParams americanCall = kCall;
    americanCall.exerciseType = ExerciseType::AMERICAN;
    checkNear(cn.price(americanCall).price, cn.price(kCall).price, 1e-6,
              "FD American call equals European call (q=0)");

    // Deep in-the-money American put must be worth at least intrinsic.
    OptionParams deepPut(50.0, 100.0, 0.05, 0.20, 1.0, 0.0, OptionType::PUT,
                         ExerciseType::AMERICAN);
    check(cn.price(deepPut).price >= 50.0 - 1e-9,
          "FD deep ITM American put >= intrinsic");

    checkThrows([] { FiniteDifference(2, 100); }, "too-coarse spot grid rejected");
    checkThrows([] { FiniteDifference(100, 0); }, "zero time steps rejected");
}

void testImpliedVolatility() {
    using namespace Options;
    ImpliedVolatility iv;
    BlackScholes bs;

    // Round trip: price at a known sigma, recover that sigma. Sweep vol
    // levels, moneyness, both option types, and a dividend yield.
    const double sigmas[] = {0.05, 0.20, 0.35, 0.80};
    for (const double sigma : sigmas) {
        OptionParams p = kCall;
        p.volatility = sigma;
        const double price = bs.price(p).price;
        checkNear(iv.solve(price, p).impliedVol, sigma, 1e-6,
                  "IV round trip call sigma=" + std::to_string(sigma));
    }
    {
        OptionParams p(80.0, 100.0, 0.05, 0.25, 0.5, 0.03, OptionType::PUT,
                       ExerciseType::EUROPEAN);
        const double price = bs.price(p).price;
        checkNear(iv.solve(price, p).impliedVol, 0.25, 1e-6,
                  "IV round trip OTM-spot put with dividends");
    }

    // The reference option's textbook price must invert to its 20% vol.
    checkNear(iv.solve(kBsCall, kCall).impliedVol, 0.20, 1e-5,
              "IV inverts reference call price");
    checkNear(iv.solve(kBsPut, kPut).impliedVol, 0.20, 1e-5,
              "IV inverts reference put price");

    // Deep in-the-money, short maturity: vega is tiny and Newton has almost
    // nothing to work with. The solver must still recover the vol (via the
    // Brent fallback when needed).
    {
        OptionParams p(100.0, 60.0, 0.05, 0.30, 0.25, 0.0, OptionType::CALL,
                       ExerciseType::EUROPEAN);
        const double price = bs.price(p).price;
        checkNear(iv.solve(price, p).impliedVol, 0.30, 1e-4,
                  "IV low-vega deep ITM call");
    }

    // Prices outside the no-arbitrage band have no implied vol.
    checkThrows([&] { (void)iv.solve(2.0, kCall); },
                "IV rejects price below intrinsic bound");
    checkThrows([&] { (void)iv.solve(100.0, kCall); },
                "IV rejects price above spot bound");
    checkThrows(
        [&] {
            OptionParams american = kCall;
            american.exerciseType = ExerciseType::AMERICAN;
            (void)iv.solve(10.0, american);
        },
        "IV rejects American options");
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
    testLongstaffSchwartz();
    testAsianOptions();
    testBarrierOptions();
    testFiniteDifference();
    testImpliedVolatility();
    testConvergenceAnalyzer();
    testPerformanceBenchmark();

    std::cout << (checks - failures) << "/" << checks << " checks passed\n";
    return failures;
}
