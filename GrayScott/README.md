# Advection-Diffusion-Reaction (Gray-Scott) Example

$$\frac{\partial y}{\partial t} + \vec{a} \cdot \nabla y -  \nabla \cdot ( D \nabla y ) = R$$

where $$y = [u v]$$ is the concentration of chemical species $$u$$ and $$v$$,
$$\vec{a}$$ is the advection speed, $$D$$ is the diffusion coefficient, and R is
the reaction term:

$$R_u = - u v^2 + A (1-u)$$

$$R_v =   u v^2 - (A + B) v$$

## Problem Options

| Option          | Type   | Description                                        | Default  |
| ----------------|--------|----------------------------------------------------|----------|
| `n_cell`        | `int`  | number of cells on each side of the square domain  | 256      |
| `max_grid_size` | `int`  | max size of boxes in box array                     | 64       |
| `plot_int`      | `int`  | enable (1) or disable (0) plots                    | 0        |
| `stepper`       | `int`  | use CVODE (0), ARKStep (1), MRIStep (2)            | 0        |
| `cvode_method`  | `int`  | use BDF (0) or Adams (1) methods in CVODE          | 0        |
| `arkode_order`  | `int`  | ARKStep method order                               | 4        |
| `nls_method`    | `int`  | use Newton (0) or fixed-point (1) solver           | 0        | 
| `nls_max_iter`  | `int`  | maximum number of nonlinear iterations             | 3        | 
| `nls_fp_acc`    | `int`  | number of fixed-point acceleration vectors         | 0        | 
| `ls_max_iter`   | `int`  | maximum number of linear iterations                | 5        | 
| `rtol`          | `Real` | relative tolerance                                 | 1e-4     | 
| `atol`          | `Real` | absolute tolerance                                 | 1e-9     | 
| `tfinal`        | `Real` | final integration time                             | 1e4      | 
| `dtout`         | `Real` | output frequency                                   | `tfinal` | 
| `advCoeffx`     | `Real` | advection speed in the x-direction                 | 5e-4     |
| `advCoeffy`     | `Real` | advection speed in the y-direction                 | 5e-4     |
| `diffCoeff`     | `Real` | diffusion coefficient                              | 2e-5     |
