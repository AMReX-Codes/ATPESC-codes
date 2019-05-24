# soruce files for GrayScott example
add_executable(Advection
  Advection.cpp
  $<TARGET_OBJECTS:nvector_multifab>
  $<TARGET_OBJECTS:diffop>)

# shared include directory
target_include_directories(Advection
  PRIVATE
  ${SHARED_INCLUDE_DIR})

# link to external libraries
target_link_libraries(Advection
  AMReX::amrex
  SUNDIALS::arkode
  SUNDIALS::cvode
  ${EXTRA_LINK_LIBS})

