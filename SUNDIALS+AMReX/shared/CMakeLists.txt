# object library for NVector wrapper
add_library(nvector_multifab OBJECT NVector_Multifab.cpp)

target_link_libraries(nvector_multifab
  PRIVATE
  AMReX::amrex
  SUNDIALS::cvode)

# object library for differential operators
add_library(diffop OBJECT DiffOp2D.cpp)

target_link_libraries(diffop
  PRIVATE
  AMReX::amrex)

# object library for reactions
add_library(reactions OBJECT Reactions.cpp)

target_link_libraries(reactions
  PRIVATE
  AMReX::amrex)

# path to shared include files
set(SHARED_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
