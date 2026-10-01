cmake_minimum_required(VERSION 3.23)

# macOS: the new LAPACK of Accelerate needs macOS 13.3 or newer. This has to be
# set before project(), it is ignored on other systems.
set(CMAKE_OSX_DEPLOYMENT_TARGET "13.3" CACHE STRING "Minimum macOS version")

project(Quanty LANGUAGES CXX)

# Optimise by default (-O3). Do not add -ffast-math or -Ofast, they cost accuracy.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Use C++11 as requested by the project
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compilation options
option(USE_INTEL "Use Intel ICX/ICPX compiler" OFF)
option(USE_OPENMP "Enable OpenMP parallelism" ON)
option(USE_MKL "Enable Intel MKL" OFF)

if(USE_INTEL)
    find_program(ICPX_PATH icpx)
    if(ICPX_PATH)
        message(STATUS "Intel compiler (ICPX) found: ${ICPX_PATH}")
        set(CMAKE_C_COMPILER icx)
        set(CMAKE_CXX_COMPILER icpx)
    else()
        message(FATAL_ERROR "Intel compiler requested, but not found.")
    endif()
endif()

# Define sources
set(SOURCES LuaInterface/main.cpp)

# Define the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Ensure the use of MINGW/MSYS on Windows
if(WIN32 AND NOT MINGW)
  message(FATAL_ERROR "MSVC is not supported for this project. Please use MinGW/MSYS on Windows.")
endif()

# Detect platform and define appropriate LUA_USE_XXX macro
if(WIN32)
  message(STATUS "Detected platform: Windows")
elseif(APPLE)
  message(STATUS "Detected platform: macOS")
  target_compile_definitions(${PROJECT_NAME} PRIVATE LUA_USE_MACOSX)
  # LUA_USE_MACOSX turns on the readline support of the Lua prompt
  target_link_libraries(${PROJECT_NAME} PRIVATE readline)
elseif(UNIX)
  message(STATUS "Detected platform: Linux/Unix")
  target_compile_definitions(${PROJECT_NAME} PRIVATE LUA_USE_LINUX)
endif()

# Optionally enable OpenMP
if(USE_OPENMP)
  # macOS: Apple's clang has no OpenMP library of its own. Use the static libomp
  # of Homebrew (brew install libomp), such that Quanty also runs on Macs
  # without Homebrew. Set OpenMP_libomp_LIBRARY and OpenMP_CXX_INCLUDE_DIR to use
  # another libomp, for example a universal one.
  if(APPLE)
    find_file(OpenMP_libomp_LIBRARY libomp.a
              PATHS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
    find_path(OpenMP_CXX_INCLUDE_DIR omp.h
              PATHS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include)
  endif()
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
  else()
    message(WARNING "OpenMP not found or not supported by the compiler")
  endif()
endif()

# Optional MKL support
if(USE_MKL)
  find_package(MKL CONFIG)
  if(MKL_FOUND)
    message(STATUS "MKL found and will be used")
    target_link_libraries(${PROJECT_NAME} PRIVATE mkl_intel_lp64 mkl_sequential mkl_core)
  else()
    message(WARNING "MKL requested but not found. Falling back to BLAS/LAPACK.")
  endif()
endif()

# Fallback to standard BLAS/LAPACK if MKL is not found or not requested
if(NOT USE_MKL OR NOT MKL_FOUND)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  if(BLAS_FOUND AND LAPACK_FOUND)
    message(STATUS "Using BLAS: ${BLAS_LIBRARIES}")
    message(STATUS "Using LAPACK: ${LAPACK_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} m pthread dl)
    # macOS: the legacy LAPACK of Accelerate has bugs (TestScripts/Matrix/LapackBug.Quanty),
    # use its new LAPACK instead
    if(APPLE AND LAPACK_LIBRARIES MATCHES "Accelerate")
      target_compile_definitions(${PROJECT_NAME} PRIVATE ACCELERATE_NEW_LAPACK)
    endif()
  else()
    message(FATAL_ERROR "Could not find BLAS or LAPACK")
  endif()
endif()


# Enable Fortran BLAS routines
target_compile_definitions(${PROJECT_NAME} PRIVATE BLASunderscore)

# Link to system libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES}
                                              ${BLAS_LIBRARIES} m pthread dl)
