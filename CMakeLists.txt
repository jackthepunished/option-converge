cmake_minimum_required(VERSION 3.15)
project(OptionsConvergenceEngine VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags for maximum performance
if(MSVC)
    # MSVC optimization flags
    add_compile_options(
        /W4                    # Warning level 4
        /O2                    # Maximize speed
        /Oi                    # Enable intrinsic functions
        /Ot                    # Favor fast code
        /GL                    # Whole program optimization
        /fp:fast               # Fast floating point model
        /arch:AVX2             # Enable AVX2 instructions (if supported)
    )
    # Linker flags for MSVC
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LTCG")
    set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} /LTCG")
else()
    # GCC/Clang optimization flags
    add_compile_options(
        -Wall                  # Enable all warnings
        -Wextra                # Extra warnings
        -Wpedantic             # Pedantic warnings
        -O3                    # Maximum optimization
        -march=native          # Optimize for current CPU
        -mtune=native          # Tune for current CPU
        -ffast-math            # Fast math operations
        -funroll-loops         # Unroll loops
        -finline-functions     # Inline functions
        -flto                  # Link-time optimization
    )
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
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
    # src/BinomialTree.cpp
    # src/MonteCarlo.cpp
    # src/ConvergenceAnalyzer.cpp
    # src/PerformanceBenchmark.cpp
)

# Create library
add_library(optionslib ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(optionslib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Main executable (commented out until main.cpp is implemented)
# add_executable(option_pricer src/main.cpp)
# target_link_libraries(option_pricer optionslib)

# Test executable (commented out until tests are implemented)
# add_executable(option_tests tests/test_main.cpp)
# target_link_libraries(option_tests optionslib)

# Enable testing
# enable_testing()
# add_test(NAME OptionTests COMMAND option_tests)

# Installation
# install(TARGETS option_pricer DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
