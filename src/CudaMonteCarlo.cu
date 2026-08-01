#include "options/CudaMonteCarlo.h"

#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace Options {

namespace {

constexpr int kBlockSize = 256;

void cudaCheck(cudaError_t err, const char* what) {
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA error in ") + what + ": " +
                                 cudaGetErrorString(err));
    }
}

// One thread, one path. Philox is counter-based: seeding with the path index
// as the subsequence gives every path its own stream regardless of how the
// grid is shaped, which is what makes the estimate launch-geometry
// independent. Block-level shared-memory reduction, then one atomicAdd per
// block into the global payoff moments.
__global__ void simulatePaths(double spot, double strike, double drift, double vol,
                              double discount, bool isCall, unsigned long long seed,
                              size_t numPaths, size_t numSteps,
                              double* sum, double* sumSq) {
    const size_t path = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;

    double payoff = 0.0;
    if (path < numPaths) {
        curandStatePhilox4_32_10_t state;
        curand_init(seed, path, 0, &state);

        double logS = log(spot);
        for (size_t step = 0; step < numSteps; ++step) {
            logS += drift + vol * curand_normal_double(&state);
        }
        const double sT = exp(logS);
        const double intrinsic = isCall ? sT - strike : strike - sT;
        payoff = discount * fmax(intrinsic, 0.0);
    }

    __shared__ double blockSum[kBlockSize];
    __shared__ double blockSumSq[kBlockSize];
    blockSum[threadIdx.x] = payoff;
    blockSumSq[threadIdx.x] = payoff * payoff;
    __syncthreads();

    for (int stride = kBlockSize / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < static_cast<unsigned>(stride)) {
            blockSum[threadIdx.x] += blockSum[threadIdx.x + stride];
            blockSumSq[threadIdx.x] += blockSumSq[threadIdx.x + stride];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        atomicAdd(sum, blockSum[0]);
        atomicAdd(sumSq, blockSumSq[0]);
    }
}

} // namespace

CudaMonteCarlo::CudaMonteCarlo(size_t numPaths, size_t numSteps, unsigned long long seed)
    : numPaths_(numPaths), numSteps_(numSteps), seed_(seed) {
    if (numPaths_ == 0) throw std::invalid_argument("CUDA MC needs at least one path");
    if (numSteps_ == 0) throw std::invalid_argument("CUDA MC needs at least one step");
}

bool CudaMonteCarlo::available() noexcept {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

CudaMonteCarlo::Result CudaMonteCarlo::price(const OptionParams& params) const {
    if (!params.isEuropean()) {
        throw std::invalid_argument("CUDA Monte Carlo prices European options only.");
    }
    if (params.hasDiscreteDividends()) {
        throw std::invalid_argument("CUDA Monte Carlo does not model discrete dividends.");
    }
    if (!available()) {
        throw std::runtime_error("No usable CUDA device.");
    }

    const double dt = params.timeToMaturity / static_cast<double>(numSteps_);
    const double drift =
        (params.riskFreeRate - params.dividendYield - 0.5 * params.volatility * params.volatility) *
        dt;
    const double vol = params.volatility * std::sqrt(dt);
    const double discount = std::exp(-params.riskFreeRate * params.timeToMaturity);

    double* device = nullptr;
    cudaCheck(cudaMalloc(&device, 2 * sizeof(double)), "cudaMalloc");
    cudaCheck(cudaMemset(device, 0, 2 * sizeof(double)), "cudaMemset");

    const auto blocks =
        static_cast<unsigned>((numPaths_ + kBlockSize - 1) / kBlockSize);
    simulatePaths<<<blocks, kBlockSize>>>(params.spotPrice, params.strikePrice, drift, vol,
                                          discount, params.isCall(), seed_, numPaths_,
                                          numSteps_, device, device + 1);

    double host[2] = {0.0, 0.0};
    const cudaError_t copyErr =
        cudaMemcpy(host, device, sizeof host, cudaMemcpyDeviceToHost);
    cudaFree(device);
    cudaCheck(copyErr, "cudaMemcpy");
    cudaCheck(cudaGetLastError(), "kernel launch");

    const double n = static_cast<double>(numPaths_);
    const double mean = host[0] / n;
    const double variance = std::fmax(host[1] / n - mean * mean, 0.0);
    return {mean, std::sqrt(variance / n)};
}

} // namespace Options
