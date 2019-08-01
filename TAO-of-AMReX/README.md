modified from amrex/Tutorials/LinearSolvers/ABecLaplacian_C example

- note: P, G have the same AMReX-defined ordering, number of entries for BC values on each process
- note: bcs: dirichlet on top, bottom, left. Neumann=0 on right. Solution to target is right edge w/o corners.

TODO:
- Fix ext BC fill indexing.

