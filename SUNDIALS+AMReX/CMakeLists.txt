cmake_minimum_required(VERSION 3.17)
cmake_policy(SET CMP0109 NEW)

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
option(ENABLE_HIP "Enable HIP" OFF)

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

if(ENABLE_HIP)
  if(NOT DEFINED ROCM_PATH)
    if(NOT DEFINED ENV{ROCM_PATH})
      set(ROCM_PATH "/opt/rocm/" CACHE PATH "Path to which ROCm has been installed")
    else()
      set(ROCM_PATH "$ENV{ROCM_PATH}" CACHE PATH "Path to which ROCm has been installed")
    endif()
  endif()

  if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
      set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
    else()
      set(HIP_PATH "$ENV{HIP_PATH}" CACHE PATH "Path to which HIP has been installed")
    endif()
  endif()

  if(NOT DEFINED HIP_PLATFORM)
    if(NOT DEFINED ENV{HIP_PLATFORM})
      set(HIP_PLATFORM "hcc" CACHE STRING "HIP platform (hcc, nvcc)")
    else()
      set(HIP_PLATFORM "$ENV{HIP_PLATFORM}" CACHE STRING "HIP platform (hcc, nvcc)")
    endif()
  endif()

  # Set CMAKE_PREFIX_PATH as the hip-config.cmake has some find_package calls
  # which don't have the proper path set (not sure if this is a bug or
  # intentional), so without this they will fail even if we provide the PATH
  # option to find_package(HIP).
  set(CMAKE_PREFIX_PATH "${ROCM_PATH};${HIP_PATH}")
  find_package(HIP REQUIRED CONFIG)
  find_package(rocrand REQUIRED CONFIG)
  find_package(rocprim REQUIRED CONFIG)
  find_package(hiprand REQUIRED CONFIG)

  if("${HIP_COMPILER}" STREQUAL "hcc")
    print_error("Deprecated HCC compiler is not supported" "Please update ROCm")
  endif()

  message(STATUS "HIP version:      ${HIP_VERSION}")
  message(STATUS "HIP platform:     ${HIP_PLATFORM}")
  message(STATUS "HIP compiler:     ${HIP_COMPILER}")
  message(STATUS "HIP linker:       ${CMAKE_CXX_LINK_EXECUTABLE}")
  message(STATUS "AMD targets:      ${AMDGPU_TARGETS}")
endif()

# Check for AMReX

# Determine AMReX components needed
if(NOT AMREX_FIND_COMPONENTS)
  set(AMREX_FIND_COMPONENTS "SUNDIALS")
  if(ENABLE_CUDA)
    list(APPEND AMREX_FIND_COMPONENTS "CUDA")
  endif()
  if(ENABLE_HIP)
    list(APPEND AMREX_FIND_COMPONENTS "HIP")
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
  if(ENABLE_HIP)
    list(APPEND SUNDIALS_FIND_COMPONENTS "nvechip")
  endif()
endif()

find_package(SUNDIALS REQUIRED COMPONENTS ${SUNDIALS_FIND_COMPONENTS}
  HINTS ${SUNDIALS_DIR} ${SUNDIALS_ROOT} ${CMAKE_PREFIX_PATH}
  NO_DEFAULT_PATH)

# Check for math library
find_library(MATH_LIBRARY NAMES m)

# Targets to build
set(targets
  "HandsOn1"
  "HandsOn2"
  "HandsOn3")

# Add the build and install targets for each example
foreach(target ${targets})

  if(ENABLE_CUDA)
    set(target_name "${target}.CUDA.exe")
  elseif(ENABLE_HIP)
    set(target_name "${target}.HIP.exe")
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
    AMReX::amrex_2d
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
  install(TARGETS ${target_name} DESTINATION ${CMAKE_INSTALL_PREFIX})

endforeach()

# Create a symlink in the build directory to each input file in source directory
file(GLOB FILES_TO_LINK "${CMAKE_CURRENT_SOURCE_DIR}/inputs*")
foreach(src ${FILES_TO_LINK})
  get_filename_component(filename ${src} NAME)
  set(link_location "${CMAKE_CURRENT_BINARY_DIR}/${filename}")
  if(EXISTS ${link_location})
    file(REMOVE ${link_location})
  endif()
  file(CREATE_LINK "${src}" "${link_location}" SYMBOLIC)
endforeach()

# Install input files and postprocessing scripts
install(FILES
  diag_tools.py
  inputs
  inputs-1
  inputs-2
  inputs-3
  inputs-ref
  DESTINATION ${CMAKE_INSTALL_PREFIX})

install(FILES
  process_ARKStep_diags.py
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  DESTINATION ${CMAKE_INSTALL_PREFIX})

# Install reference solution
install(DIRECTORY
  reference_solution
  DESTINATION ${CMAKE_INSTALL_PREFIX})

# find_program(AMREX_FCOMPARE fcompare ${AMREX_ROOT}/bin)
# if(AMREX_FCOMPARE)
#   message(STATUS "Found fcompare: ${AMREX_FCOMPARE}")
#   install(FILES
#     ${AMREX_FCOMPARE}
#     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
#     DESTINATION ${CMAKE_INSTALL_PREFIX})
# else()
#   message(STATUS "Not found fcompare!")
# endif()

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
endif()
