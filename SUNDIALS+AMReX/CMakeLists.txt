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

# Determine AMReX components needed
if(NOT AMREX_FIND_COMPONENTS)
  set(AMREX_FIND_COMPONENTS "SUNDIALS")
  if(ENABLE_CUDA)
    list(APPEND AMREX_FIND_COMPONENTS "CUDA")
  endif()
endif()

find_package(AMReX REQUIRED COMPONENTS ${AMREX_FIND_COMPONENTS}
  HINTS ${AMREX_DIR} ${AMREX_ROOT} ${CMAKE_PREFIX_PATH}
  NO_DEFAULT_PATH)

# Check for SUNDIALS

# Determine SUNDIALS components needed
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

# Targets to build
set(targets
  "HandsOn1"
  "HandsOn2"
  "HandsOn3")

# Add the build and install targets for each example
foreach(target ${targets})

  if(ENABLE_CUDA)
    set(target_name "${target}.CUDA.exe")
  else()
    set(target_name "${target}.exe")
  endif()

  add_executable(${target_name} ${target}.cpp HandsOn_main.cpp)

  if(ENABLE_CUDA)
    set_source_files_properties(${target}.cpp HandsOn_main.cpp
      PROPERTIES LANGUAGE CUDA)
    # Add -dc flag
    set_target_properties(${target_name}
      PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

  # Link to required libraries
  target_link_libraries(${target_name} PRIVATE
    AMReX::amrex
    SUNDIALS::arkode
    SUNDIALS::cvode
    MPI::MPI_CXX
    m)

  # Link to HYPRE using target defined by AMReX
  if (AMReX_HYPRE)
    target_link_libraries(${target_name} PRIVATE HYPRE)
  endif()

  # Link to libraries needed for CUDA support
  if(ENABLE_CUDA)
    target_link_libraries(${target_name} PRIVATE
      SUNDIALS::nveccuda
      CUDA::cusparse
      CUDA::curand
      CUDA::cublas
      CUDA::nvToolsExt)
  endif()

  # Install execuatable
  install(TARGETS ${target_name} RUNTIME DESTINATION bin)

endforeach()

# Install input files and postprocessing scripts
install(FILES
  diag_tools.py
  inputs
  inputs-1
  inputs-2
  inputs-3
  inputs-ref
  DESTINATION bin)

install(FILES
  process_ARKStep_diags.py
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  DESTINATION bin)

# Install reference solution
install(DIRECTORY
  reference_solution
  DESTINATION bin)

# Add target to clean up plots when running in the build directory
add_custom_target(pltclean
  COMMAND ${CMAKE_COMMAND} -E rm -rf plt*/ *.png *_diagnostics.txt __pycache__)

# Add target with list of plot headers for making a movie in Visit
if(UNIX)
  add_custom_target(movie
    COMMAND ls -1 plt*/Header | tee movie.visit)
endif()

# Install source files
if(INSTALL_SOURCE)
  install(FILES
    CMakeLists.txt
    HandsOn.hpp
    HandsOn1.cpp
    HandsOn2.cpp
    HandsOn3.cpp
    HandsOn_main.cpp
    README.md
    diag_tools.py
    inputs
    inputs-1
    inputs-2
    inputs-3
    inputs-ref
    DESTINATION SUNDIALS+AMReX)

  install(FILES
    process_ARKStep_diags.py
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    DESTINATION bin)
endif()
