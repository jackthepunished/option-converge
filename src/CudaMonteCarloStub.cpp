// Compiled instead of CudaMonteCarlo.cu when CMake finds no CUDA toolchain.
// Keeps the class linkable everywhere; pricing honestly refuses.

#include "options/CudaMonteCarlo.h"

#include <stdexcept>

namespace Options {

CudaMonteCarlo::CudaMonteCarlo(size_t numPaths, size_t numSteps, unsigned long long seed)
    : numPaths_(numPaths), numSteps_(numSteps), seed_(seed) {
    if (numPaths_ == 0) throw std::invalid_argument("CUDA MC needs at least one path");
    if (numSteps_ == 0) throw std::invalid_argument("CUDA MC needs at least one step");
}

bool CudaMonteCarlo::available() noexcept {
    return false;
}

CudaMonteCarlo::Result CudaMonteCarlo::price(const OptionParams&) const {
    throw std::runtime_error(
        "This build has no CUDA support; rebuild with a CUDA toolchain to price on the GPU.");
}

} // namespace Options
