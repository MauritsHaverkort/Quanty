# Building Quanty

Quanty is intended to compile on all Linux, Windows and macOS systems found in the wild. This includes a ten-year-old cluster at some university institute that runs fine and has never been updated. We cannot, and do not, assume that libraries are available. Quanty links to BLAS and LAPACK, but should not depend on other external libraries (i.e. no Boost, HDF5 or FFT libraries).

The code compiles on Windows, Linux and macOS with the Intel, GNU, Clang and other compilers. We follow the C++11 standard to stay compatible with different systems.

The I/O functions contain code specific to each operating system. On Linux and macOS you select it on the compile command line with

```
-D LUA_USE_LINUX
-D LUA_USE_MACOSX
```

respectively. On Windows it is selected automatically.

You only need to compile the file `LuaInterface/main.cpp`; it includes all headers and other files needed.

## Linux, Intel compiler

A possible compile command for the Intel compiler on Linux is

```bash
icc LuaInterface/main.cpp -D LUA_USE_LINUX -O3 -ipo -static -qopenmp -I"${MKLROOT}/include" -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -o Quanty
```

i.e. we link to OpenMP, BLAS, LAPACK and the standard libraries.

## Linux, GNU g++ or clang++

A possible command for GNU g++ or clang++ is

```bash
g++ LuaInterface/main.cpp -D LUA_USE_LINUX -DBLASunderscore -O3 -fopenmp -o Quanty -lblas -llapack
```

Note the `BLASunderscore` macro: Quanty uses it to support the different naming conventions of BLAS and LAPACK libraries.

## macOS

On macOS, BLAS and LAPACK come from Apple's Accelerate framework. Accelerate contains two LAPACK versions. The old one (LAPACK 3.2.1, the default) has bugs; for example, `zheev` returns wrong eigenvectors for a Hermitian matrix of dimension 130 or more whose last rows and columns are zero (see `TestScripts/Matrix/LapackBug.Quanty`). Compile with `-DACCELERATE_NEW_LAPACK` to use the new LAPACK instead. This requires macOS 13.3 or newer. Do not use `ACCELERATE_LAPACK_ILP64`; Quanty passes 32-bit integers to BLAS and LAPACK.

Apple's clang does not come with OpenMP. Install libomp (for example with Homebrew: `brew install libomp`) and link the static library, so that the binary does not depend on Homebrew on the computer it runs on. The Homebrew path is `/usr/local/opt/libomp` on Intel Macs and `/opt/homebrew/opt/libomp` on Apple silicon.

A possible command on macOS is

```bash
clang++ LuaInterface/main.cpp -D LUA_USE_MACOSX -DBLASunderscore -DACCELERATE_NEW_LAPACK -O3 -mmacosx-version-min=13.3 -Xclang -fopenmp -I/usr/local/opt/libomp/include -o Quanty -framework Accelerate -lreadline /usr/local/opt/libomp/lib/libomp.a
```

A universal binary, which runs on both Intel and Apple silicon Macs, is made by running the command above twice, once with `-arch x86_64 -o Quanty-x86_64` and once with `-arch arm64 -o Quanty-arm64` (each needs a `libomp.a` for that architecture), and joining the two with

```bash
lipo -create Quanty-x86_64 Quanty-arm64 -output Quanty
```

## CMake

`CMakeLists.txt` does all of the above on macOS by itself (new LAPACK, static libomp from Homebrew, readline, `-O3`, macOS 13.3):

```bash
cmake -S . -B build
cmake --build build
```

For a universal binary, add `-DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" -DOpenMP_libomp_LIBRARY=<path to a universal libomp.a>` to the first command.

## Compiler flags to avoid

- Do not use `-ffast-math` or `-Ofast`. They change the order of summations, flush small numbers to zero and break the tests for NaN, and thereby the numerical accuracy.
- Do not use `-march=native` for a binary that should run on other computers.
- `-flto` gains nothing, as Quanty is compiled as a single file.

## Running

The number of threads is set with the environment variable `OMP_NUM_THREADS`.
