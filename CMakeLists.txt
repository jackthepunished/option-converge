cmake_minimum_required(VERSION 3.15)
project(OptionsConvergenceEngine VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Warnings always; optimization flags only in Release so Debug builds stay
# debuggable and do not collide with MSVC's /RTC1 runtime checks.
# Fast-math relaxes IEEE semantics (reassociation, no NaN/inf assumptions).
# Overridable so the numerical impact can be measured against strict FP;
# see the README's Performance section for the measured trade-off.
option(OPTIONS_FAST_MATH "Enable -ffast-math / /fp:fast in Release builds" ON)

if(MSVC)
    add_compile_options(
        /W4                    # Warning level 4
        "$<$<CONFIG:Release>:/O2;/Oi;/Ot;/GL;/arch:AVX2>"
    )
    if(OPTIONS_FAST_MATH)
        add_compile_options("$<$<CONFIG:Release>:/fp:fast>")
    endif()
    # /LTCG pairs with /GL, so Release only
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} /LTCG")
else()
    add_compile_options(
        -Wall                  # Enable all warnings
        -Wextra                # Extra warnings
        -Wpedantic             # Pedantic warnings
        "$<$<CONFIG:Release>:-O3;-funroll-loops;-finline-functions>"
    )
    if(OPTIONS_FAST_MATH)
        add_compile_options("$<$<CONFIG:Release>:-ffast-math>")
    endif()
    # Not every compiler accepts -march=native (e.g. older Apple Clang on arm64)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        add_compile_options("$<$<CONFIG:Release>:-march=native;-mtune=native>")
    endif()
    # LTO only where the entire toolchain (compiler, archiver, linker plugin)
    # supports it. Hardcoding -flto breaks when the compiler and the linker's
    # LTO plugin come from different versions, as on CI runners where clang 18
    # emits bitcode the system bfd/LLVMgold plugin (LLVM 16) cannot read.
    include(CheckIPOSupported)
    check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_CHECK_OUTPUT LANGUAGES CXX)
    if(LTO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
        message(STATUS "IPO/LTO enabled for Release")
    else()
        message(STATUS "IPO/LTO not supported by this toolchain, continuing without it")
    endif()
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Find OpenMP for parallel processing (optional but recommended)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found - parallel processing enabled")
endif()

# Source files
set(SOURCES
    src/Option.cpp
    src/PricingEngine.cpp
    src/BlackScholes.cpp
    src/AsianOption.cpp
    src/BarrierOption.cpp
    src/BinomialTree.cpp
    src/Calibration.cpp
    src/DigitalOption.cpp
    src/FiniteDifference.cpp
    src/Heston.cpp
    src/ImpliedVolatility.cpp
    src/LongstaffSchwartz.cpp
    src/LookbackOption.cpp
    src/MonteCarlo.cpp
    src/ConvergenceAnalyzer.cpp
    src/PerformanceBenchmark.cpp
)

# Optional CUDA support: build the GPU Monte Carlo engine when a CUDA
# toolchain exists, otherwise compile a stub that reports unavailability.
# CI runners without nvcc take the stub path and stay green.
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
            set(CMAKE_CUDA_ARCHITECTURES native)
        else()
            set(CMAKE_CUDA_ARCHITECTURES 75)
        endif()
    endif()
    list(APPEND SOURCES src/CudaMonteCarlo.cu)
    message(STATUS "CUDA toolchain found - GPU Monte Carlo enabled")
else()
    list(APPEND SOURCES src/CudaMonteCarloStub.cpp)
    message(STATUS "No CUDA toolchain - GPU Monte Carlo built as stub")
endif()

# Create library
add_library(optionslib ${SOURCES})
if(CMAKE_CUDA_COMPILER)
    target_compile_definitions(optionslib PUBLIC OPTIONS_HAS_CUDA)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(optionslib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Main executable
add_executable(option_pricer src/main.cpp)
target_link_libraries(option_pricer optionslib)

# Test executable
add_executable(option_tests tests/test_main.cpp)
target_link_libraries(option_tests optionslib)

# Enable testing
enable_testing()
add_test(NAME OptionTests COMMAND option_tests)

# Installation
# install(TARGETS option_pricer DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
