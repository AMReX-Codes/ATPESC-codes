cmake_minimum_required(VERSION 3.17)

# Prevent in-source builds
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(FATAL_ERROR "\nIn-source builds are prohbited!\n")
endif ()

# ATPESC example project
project(atpesc-sundials-amrex-example
  DESCRIPTION "ATPESC SUNDIALS+AMReX Example"
  LANGUAGES C CXX)

# Build options
option(ENABLE_CUDA "Enable CUDA" OFF)

set(AMREX_ROOT "$ENV{AMREX_ROOT}"
  CACHE PATH "Path to AMReX installation directory")

set(SUNDIALS_ROOT "$ENV{SUNDIALS_DIR}"
  CACHE PATH "Path to SUNDIALS installation directory")

# Check for MPI
find_package(MPI REQUIRED)

# Compiler options
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
endif()

# Setup CUDA
if(ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  find_package(CUDAToolkit REQUIRED)
endif()

# Check for AMReX

# determine AMReX components needed
if(NOT AMREX_FIND_COMPONENTS)
  set(AMREX_FIND_COMPONENTS )
  if(ENABLE_CUDA)
    list(APPEND AMREX_FIND_COMPONENTS "CUDA")
  endif()
endif()

find_package(AMReX REQUIRED COMPONENTS ${AMREX_FIND_COMPONENTS}
  HINTS ${AMREX_DIR} ${AMREX_ROOT} ${CMAKE_PREFIX_PATH}
  NO_DEFAULT_PATH)

# Check for SUNDIALS

# determine SUNDIALS components needed
if(NOT SUNDIALS_FIND_COMPONENTS)
  set(SUNDIALS_FIND_COMPONENTS
    "arkode"
    "cvode")
  if(ENABLE_CUDA)
    list(APPEND SUNDIALS_FIND_COMPONENTS "nveccuda")
  endif()
endif()

find_package(SUNDIALS REQUIRED COMPONENTS ${SUNDIALS_FIND_COMPONENTS}
  HINTS ${SUNDIALS_DIR} ${SUNDIALS_ROOT} ${CMAKE_PREFIX_PATH}
  NO_DEFAULT_PATH)

# Check for math library
find_library(MATH_LIBRARY NAMES m)

# Don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# When building, don't use the install RPATH already (but later on when
# installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")

# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif()

# Targets
set(targets
  "HandsOn1"
  "HandsOn2"
  "HandsOn3")

# Add the build and install targets for each example
foreach(target ${targets})

  add_executable(${target}.exe ${target}.cpp HandsOn_main.cpp)

  if(ENABLE_CUDA)
    set_source_files_properties(${target}.cpp HandsOn_main.cpp
      PROPERTIES LANGUAGE CUDA)
    # Add -dc flag
    set_target_properties(${target}.exe
      PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

  target_link_libraries(${target}.exe PRIVATE
    AMReX::amrex
    SUNDIALS::arkode
    SUNDIALS::cvode
    MPI::MPI_CXX
    m)

  if(ENABLE_CUDA)
    target_link_libraries(${target}.exe PRIVATE
      SUNDIALS::nveccuda
      CUDA::cusparse
      CUDA::curand
      CUDA::nvToolsExt)
  endif()

  install(TARGETS ${target}.exe RUNTIME DESTINATION bin)

endforeach()

# Add target to clean up plots when running in the build directory
add_custom_target(pltclean
  COMMAND ${CMAKE_COMMAND} -E rm -rf plt*/ *.png *_diagnostics.txt __pycache__)

if (UNIX)
  add_custom_target(movie
    COMMAND ls -1 plt*/Header | tee movie.visit)
endif()
