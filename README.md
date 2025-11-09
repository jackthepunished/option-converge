Option Converge

Overview
--------
This repository contains a small C++ project for option pricing and convergence analysis. It includes implementations of pricing engines, Monte Carlo and Binomial Tree methods, utilities for benchmarking, and tools to compare numerical convergence. The project is configured with CMake and organized to separate headers, sources, tests, and result artifacts.

Key goals
- Provide clear, maintainable C++ code for option pricing methods.
- Allow reproducible benchmarking and convergence analysis.
- Keep the build process simple using CMake.

Repository layout
---------------
- `CMakeLists.txt` - Top-level CMake file used to configure and build the project.
- `include/` - Public headers. Notable subfolders:
  - `options/` - Option pricing interfaces and implementations (e.g. `BinomialTree.h`, `MonteCarlo.h`, `BlackScholes.h`).
- `src/` - Project source files (implementations).
- `tests/` - Unit and integration tests.
- `scripts/` - Helper scripts for running analyses or processing results.
- `results/` - Output directory for generated benchmarking and convergence results.
- `docs/` - Documentation and notes.

Prerequisites
-------------
- A modern C++ compiler supporting C++17 or later (for example, Clang or GCC).
- CMake (version 3.10 or newer recommended).
- Make or a similar build tool (CMake will choose an appropriate generator on macOS).

Build and run (macOS / zsh)
--------------------------
1. Create an out-of-source build directory and configure the project with CMake:

   mkdir -p build
   cmake -S . -B build

2. Build the project (parallel build recommended):

   cmake --build build -- -j$(sysctl -n hw.ncpu)

3. Run tests (if tests are configured to be discovered by CTest):

   cmake --build build --target test
   # or
   cd build && ctest --output-on-failure

4. Executables and example programs (if present) will be located in `build/` or a subfolder specified by the CMake configuration. Check `CMakeLists.txt` to confirm exact targets and installed paths.

Notes about dependencies and configuration
-----------------------------------------
- If the project adds external dependencies through `find_package` or submodules, ensure those packages are installed or available on your system prior to configuration.
- If you need to change compiler flags or enable additional features, edit `CMakeLists.txt` or pass options to `cmake` when configuring. For example:

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

Running benchmarks and generating results
----------------------------------------
- Use the binaries in `build/` or any helper scripts in `scripts/` to run benchmarks.
- Benchmark outputs and analysis results are placed in `results/`. Make sure that the `results/` directory is writable by the user running the program.

Tests
-----
- Unit tests are in the `tests/` directory. The project is expected to run tests via CTest when `make test` or the `test` target is invoked by CMake.
- If the tests rely on extra data or environment variables, those will be documented in `tests/` or the relevant source headers; check the test code for specifics.

Contributing
------------
- If you want to contribute, please fork the repository, create a topic branch, and open a pull request with a clear description of the changes.
- Keep changes small and focused. Add or update tests for new behavior.
- Follow the existing code style and conventions used across the headers in `include/` and implementation files in `src/`.

Reporting issues
----------------
- Use the repository issue tracker to report bugs or feature requests. Provide a concise description, steps to reproduce, and any relevant environment details (compiler version, OS, CMake version).

Contact
-------
For questions about the project, use the repository communication channels (issues or the project README/contact section if present) and include relevant details about your environment and what you are trying to achieve.

Acknowledgements
-----------------
This project was organized to keep experiments and numerical analysis reproducible and easy to run. If you adapt or extend the work, please mention the repository in your project notes.
