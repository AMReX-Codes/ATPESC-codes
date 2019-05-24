cmake_minimum_required(VERSION 3.10)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source build prohibited.")
endif()

project(ATPESC-examples CXX)

# check for AMREX location
if(NOT AMREX_DIR)
  set(AMREX_DIR "" CACHE FILEPATH "AMReX install directory")
  message(FATAL_ERROR "AMREX_DIR is not set")
endif()

# find AMReX
find_package(AMREX
  NAMES AMREX
  PATHS ${AMREX_DIR} ${AMREX_DIR}/lib/cmake/AMReX
  NO_DEFAULT_PATH)

if(AMREX_FOUND)
  message(STATUS "Looking for AMReX... success")
else()
  message(FATAL_ERROR "Looking for AMReX... failed")
endif()

# check for SUNDIALS location
if(NOT SUNDIALS_DIR)
  set(SUNDIALS_DIR "" CACHE FILEPATH "SUNDIALS install directory")
  message(FATAL_ERROR "SUNDIALS_DIR is not set")
endif()

# find ARKode library
find_library(ARKODE_LIBS sundials_arkode
  PATHS ${SUNDIALS_DIR}/lib ${SUNDIALS_DIR}/lib64
  NO_DEFAULT_PATH)

if(ARKODE_LIBS)
  message(STATUS "Looking for ARKode... success")
else()
  message(FATAL_ERROR "Looking for ARKode... failed")
endif()

# Create target for SUNDIALS
if(NOT TARGET SUNDIALS::arkode)
  add_library(SUNDIALS::arkode UNKNOWN IMPORTED)
  set_property(TARGET SUNDIALS::arkode PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SUNDIALS_DIR}/include)
  set_property(TARGET SUNDIALS::arkode PROPERTY IMPORTED_LOCATION ${ARKODE_LIBS})
endif()

# find CVODE library
find_library(CVODE_LIBS sundials_cvode
  PATHS ${SUNDIALS_DIR}/lib ${SUNDIALS_DIR}/lib64
  NO_DEFAULT_PATH)

if(CVODE_LIBS)
  message(STATUS "Looking for CVODE... success")
else()
  message(FATAL_ERROR "Looking for CVODE... failed")
endif()

# Create target for SUNDIALS
if(NOT TARGET SUNDIALS::cvode)
  add_library(SUNDIALS::cvode UNKNOWN IMPORTED)
  set_property(TARGET SUNDIALS::cvode PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SUNDIALS_DIR}/include)
  set_property(TARGET SUNDIALS::cvode PROPERTY IMPORTED_LOCATION ${CVODE_LIBS})
endif()

# extra link library
set(EXTRA_LINK_LIBS gfortran)

# add source directories
add_subdirectory(shared)
add_subdirectory(Advection)
add_subdirectory(GrayScott)
