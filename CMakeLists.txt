cmake_minimum_required(VERSION 3.23)
project(Quanty LANGUAGES CXX)

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
elseif(UNIX)
  message(STATUS "Detected platform: Linux/Unix")
  target_compile_definitions(${PROJECT_NAME} PRIVATE LUA_USE_LINUX)
endif()

# Optionally enable OpenMP
if(USE_OPENMP)
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
  else()
    message(FATAL_ERROR "Could not find BLAS or LAPACK")
  endif()
endif()


# Enable Fortran BLAS routines
target_compile_definitions(${PROJECT_NAME} PRIVATE BLASunderscore)

# Link to system libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES}
                                              ${BLAS_LIBRARIES} m pthread dl)
