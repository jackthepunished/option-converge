cmake_minimum_required(VERSION 3.15)
project(OptionsConvergenceEngine VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
if(MSVC)
    add_compile_options(/W4 /O2)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -O3 -march=native)
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
    src/BlackScholes.cpp
    src/BinomialTree.cpp
    src/MonteCarlo.cpp
    src/Greeks.cpp
    src/ConvergenceAnalyzer.cpp
    src/PerformanceBenchmark.cpp
    src/Utils.cpp
)

# Create library
add_library(optionslib ${SOURCES})

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
install(TARGETS option_pricer DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
