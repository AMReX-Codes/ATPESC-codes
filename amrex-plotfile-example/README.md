# AMReX Plotfile Example

This sample code initializes a MultiFab with one component "phi" using
the `(i,j,k)` indices of each cell as:

```
phi(i,j,k) = i + 100.0*j + 10000.0*k
```

The MultiFab is then written to a plotfile.

## Compiling

First, compile AMReX as a library using the directions here:
https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-libamrex

This was tested with g++ 5.4.0 with the default options supplied to `./configure` as:

```
./configure --prefix=[AMReX library prefix]
```

Then compile this example as:

```
make AMREX_LIBRARY_HOME=[AMReX library prefix]
```

## Running

If compilation is successful you should have a `main.exe` executable.

This can be run as-is with defaults corresponding to a 32^3 domain and
max grid size of 16. This will create a domain of size 32^3 cells in
3D and divide it up into 8 boxes each of size 16^3.

It will also write a plotfile named `plt_32_32_32_16`. The naming
convention for the plotfiles is:

```
plt_[# cells x]_[# cells y]_[# cellsz]_[max grid size]
```

To change these defaults, you can use command line options. For
example, this will set up a 64^3 domain and divide it up into 8 boxes
of size 32^3:

```
./main.exe n_cells=64 64 64 max_grid_size=32
```
