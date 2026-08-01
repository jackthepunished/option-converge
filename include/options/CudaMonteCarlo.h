#ifndef CUDA_MONTE_CARLO_H
#define CUDA_MONTE_CARLO_H

#include "Option.h"

#include <cstddef>

namespace Options {

// GPU Monte Carlo pricer for European vanilla options. One CUDA thread
// simulates one GBM path with a Philox counter-based generator whose
// subsequence is the path index, so the draws - and the estimate, up to
// floating-point reduction order - are independent of the launch geometry
// and reproducible for a fixed seed, matching the CPU engine's contract.
//
// The class is always declared; whether it can price depends on the build.
// When CMake finds no CUDA toolchain a stub is compiled instead, available()
// reports false, and price() throws. This keeps the library, tests, and CI
// identical in shape whether or not a GPU toolchain exists.
class CudaMonteCarlo {
public:
    struct Result {
        double price;
        double standardError;
    };

    explicit CudaMonteCarlo(size_t numPaths = 1000000, size_t numSteps = 252,
                            unsigned long long seed = 42);

    // Price a European vanilla option on the GPU. Throws for American
    // exercise, discrete dividends, or when built without CUDA.
    [[nodiscard]] Result price(const OptionParams& params) const;

    // True only when built with CUDA and at least one device is usable.
    [[nodiscard]] static bool available() noexcept;

private:
    size_t numPaths_;
    size_t numSteps_;
    unsigned long long seed_;
};

} // namespace Options

#endif // CUDA_MONTE_CARLO_H
