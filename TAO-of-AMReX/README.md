modified from amrex/Tutorials/LinearSolvers/ABecLaplacian_C example

- note: P, G have the same AMReX-defined ordering, number of entries for BC values on each process
- note: bcs: dirichlet on top, bottom, left. Neumann=0 on right. Solution to target is right edge w/o corners.

TODO:
- modify f so u_out is the right edge of the domain, not including corners
- write function to give degrees of freedom for each processor, cell ordering for boundary
- fill BC, Neumann along right edge

- Alp to give u_target function defined in physical coordinates
- Don to define u_target from the function^