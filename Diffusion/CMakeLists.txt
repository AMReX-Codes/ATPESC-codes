# soruce files for GrayScott example
add_executable(Diffusion
  Diffusion.cpp
  $<TARGET_OBJECTS:nvector_multifab>
  $<TARGET_OBJECTS:diffop>)

# shared include directory
target_include_directories(Diffusion
  PRIVATE
  ${SHARED_INCLUDE_DIR})

# link to external libraries
target_link_libraries(Diffusion
  AMReX::amrex
  SUNDIALS::arkode
  SUNDIALS::cvode
  ${EXTRA_LINK_LIBS})

# copy input file to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/inputs
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
