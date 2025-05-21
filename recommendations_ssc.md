# Suggestions for Improving the Quanty Software Infrastructure

**Date**: 07.05.2025

## Introduction

This document outlines a number of potential improvements to the Quanty code base and project infrastructure, based on our preliminary evaluation. All suggestions are meant to be non-intrusive, respectful of the existing style, and incremental, with a focus on long-term maintainability, onboarding of new contributors, and performance improvements.

## Incorporating Modern C++ Practices and [Standard Library Usage](https://en.cppreference.com/w/cpp/standard_library)

**Goal**: Improve readability, safety, and maintainability by gradually using more modern C++ features where appropriate.

- **Avoiding manual and unsafe macros**:
  - Macros like `Max()`, `Abs()`, or `Sqr()` can be replaced by standard functions ([`std::abs`](https://en.cppreference.com/w/cpp/numeric/math/abs), [`std::max`](https://en.cppreference.com/w/cpp/algorithm/max), etc.). These replacements improve safety, reduce side effects, and improve readability.
  - _Example_: `Core/BasicMath/BasicMath.cpp`
  - _Expected workload_: Low to moderate, can be done incrementally.

- **Standard math functions**:
  - Functions like `PowerInt`, `PowerDoubleInt` can be replaced with a templated function or [`std::pow`](https://en.cppreference.com/w/cpp/numeric/math/pow) where appropriate.
  - Could help avoid redundancy and errors due to inconsistent re-implementation.
  - _Example_: `Core/BasicMath/BasicMath.cpp`, e.g., `template<typename T> T power(T x, unsigned a)`
  - _Expected workload_: Low

- **Define types in headers**:
  - It is recommended to define types (like structs or classes) in header files so they can be shared across multiple source files (cpp) that can be compiled separately. This promotes code reuse and modularity, and ensures consistency and type visibility during compilation.
  - _Example_: `Core/Structdef.cpp` contains the definition of `Complex16`.
  - _Expected workload_: High

- **Pass parameters by const reference where applicable**:
  - Many helper functions pass arrays or structs (e.g., `Complex16`, fixed-length vectors) by value. This leads to unnecessary copies and makes accidental mutation harder to detect. Passing by const reference improves performance and increases safety.
  - _Example_: In `Core/BasicMath/BasicMath_SparseMatrixVector.cpp`, prefer `const SparseMatrixType& Matrix` over pass-by-value.
  - _Expected workload_: Low

- **Use of [`std::complex`](https://en.cppreference.com/w/cpp/numeric/complex)**:
  - Currently, a custom `Complex16` struct is used.
  - The STL provides a fully-featured, optimized `std::complex<double>` with arithmetic, I/O, and helper operations.
  - _Expected workload_: High, as existing code would need significant adjustment.

- **Replace manual timing with [`std::chrono`](https://en.cppreference.com/w/cpp/chrono)**:
  - The current timing system (`Core/TimeStamp.cpp`) uses manual memory allocation and timestamp management.
  - `std::chrono` is a type-safe, cross-platform, and high-resolution timing library.
  - Replace `BeginTiming("A")` / `EndTiming("A")` with a `ScopedTimer t("A")` that reports time on destruction.
  - _Expected workload_: Low

- **Remove unused functions and empty files**:
  - Unused files:
    - `Core/BasicMath/BasicMath_InterpolatingFunction.*`
    - `Core/BasicMath/BasicMath_Spline.*`
    - `Core/FiniteTemperatureLanczos.*`
    - `LuaInterface/Plot.*`
  - Cleaning these improves maintainability and reduces confusion.
  - _Expected workload_: Low

- **Replace manual I/O and string formatting with modern C++**:
  - Files like `IO/stringPrint.h` use `char*`, `sprintf`, and `realloc` for string formatting.
  - These are verbose, error-prone, and hard to maintain.
  - Use [`std::ostringstream`](https://en.cppreference.com/w/cpp/io/basic_ostringstream) or [`std::format`](https://en.cppreference.com/w/cpp/utility/format/format) (C++20), combined with [`std::string`](https://en.cppreference.com/w/cpp/string/basic_string).
  - _Expected workload_: Moderate to high

## Use of External Libraries

**Goal**: Improve performance, readability and simplify implementation of well-established numerical operations.

### Linear Algebra: [Eigen](https://eigen.tuxfamily.org)

- Use Eigen’s efficient abstractions for matrix and vector operations.
- Header-only and lightweight — no installation required.
- Optimized for modern CPUs; can delegate to BLAS/LAPACK for performance.
- Compatible with GCC, Clang, MSVC, Intel.
- Improves clarity compared to raw BLAS or custom macros.

**Comparison with alternatives**:

- **BLAS / LAPACK**: Fast, but low-level and hard to use.
- **Intel MKL**: Extremely fast but proprietary.
- **OpenBLAS / Netlib**: Open-source, but lacks abstraction.

_Example_: Replace `Dot()`, `Cross()`, `Covariance()` with `Eigen::Vector3d`, `dot()`, `cross()`.

_Expected workload_: High

### Tensors
- [Eigen](https://eigen.tuxfamily.org/dox/unsupported/eigen_tensors.html): Primarily a linear algebra library, but supports basic tensor operations through the Tensor module.
- [xtensor](https://github.com/xtensor-stack/xtensor): Multi-dimensional arrays with tensor capabilities.
- [TACO](https://github.com/tensor-compiler/taco): C++ library that computes tensor algebra expressions on sparse and dense tensors.

### Splines

- [Boost](https://live.boost.org/): Provides basic spline interpolation
- [Eigen](https://eigen.tuxfamily.org/dox/unsupported/group__Splines__Module.html): Limited support for spline interpolation
- [splinterp](https://github.com/apryor6/splinterp): fast C++ template library for parallel calculation of linear, bilinear, and trilinear interpolation from data sampled on an integer grid.

## Project Organization and Best Practices

**Goal**: Make the project easier to build, contribute to, and maintain.

### Introduce CMake-based build system

- **Cross-platform support**: Works on Linux, macOS, and Windows.
- **No dependency bloat**: Keeps build system minimal.
- **Works on legacy systems**: CMake is usually available or easily built from source.
- **Same linking behavior**: Identical compiler/linker flags can be reused.
- **Cleaner dependency management**: Centralizes defines and options.
- **IDE integration**: Enables support for VSCode, CLion, etc.
- **Improved automation**: Easier to add files, tests, and components.

_Expected workload_: Low

### Improve CI readability and maintainability

- The CI configuration currently contains many near-identical job definitions for running tests with different input scripts.
- GitLab CI supports defining job templates and using a matrix strategy to generate multiple jobs programmatically.
- This improves maintainability and reduces duplication.
- Example:

```yaml
.test-rf:
  stage: test_responsefunction
  script:
    - export OMP_NUM_THREADS=8
    - ./Quanty.exe "TestScripts/Responsefunction/${TESTFILE}.Quanty"
  dependencies:
    - build-job

test-rf-matrix:
  extends: .test-rf
  parallel:
    matrix:
      - TESTFILE: ["__eq", "__add", "__sub", "__mul"]
```

- *Effect*: Defines a reusable template and generates all jobs automatically for the listed test files.
- *Expected workload*: Low

### Coding conventions and tooling

- Adopting a unified code style (indentation, naming, etc.) aids collaboration.
- Tools like `clang-format` and pre-commit hooks could apply formatting automatically.
- These tools can be configured to match existing preferences.
- *Expected workload*: Low once initial config is defined.

### Improve project README and documentation

- README could be enhanced with a quick overview, installation instructions, usage examples.
- Especially helpful for new contributors or users who do not visit the website directly.
- Markdown is supported.
- *Expected workload*: Low
