# Advection-Diffusion Example

$$\frac{\partial u}{\partial t} + \vec{a} \cdot \nabla u -  \nabla \cdot ( D \nabla u ) = 0$$

where $$u$$ is the chemical concentration, $$\vec{a}$$ is the advection speed, and
$$D$$ is the diffusion coefficient.

## Problem Options

| Option          | Type   | Description                                        | Default  |
| ----------------|--------|----------------------------------------------------|----------|
| `n_cell`        | `int`  | number of cells on each side of the square domain  | 256      |
| `max_grid_size` | `int`  | max size of boxes in box array                     | 64       |
| `plot_int`      | `int`  | enable (1) or disable (0) plots                    | 0        |
| `stepper`       | `int`  | use CVODE (0) or ARKStep (1)                       | 0        |
| `cvode_method`  | `int`  | use BDF (0) or Adams (1) methods in CVODE          | 0        |
| `arkode_order`  | `int`  | ARKStep method order                               | 4        |
| `nls_method`    | `int`  | use Newton (0) or fixed-point (1) solver           | 0        |
| `nls_max_iter`  | `ìnt`  | maximum number of nonlinear iterations             | 3        |
| `nls_fp_acc`    | `int`  | number of fixed-point acceleration vectors         | 0        |
| `ls_max_iter`   | `int`  | maximum number of linear iterations                | 5        |
| `rhs_adv`       | `int`  | advection: disable (0), implicit (1), explicit (2) | 1        |
| `rhs_diff`      | `int`  | diffusion: disable (0), implicit (1), explicit (2) | 1        |
| `rtol`          | `Real` | relative tolerance                                 | 1e-4     |
| `atol`          | `Real` | absolute tolerance                                 | 1e-9     |
| `fixed_dt`      | `Real` | use a fixed time step size (if `fixed_dt` > 0.0)   | -1.0     |
| `tfinal`        | `Real` | final integration time                             | 1e4      |
| `dtout`         | `Real` | output frequency                                   | `tfinal` |
| `write_diag`    | `int`  | output ARKStep diagnostics to a file               | 1        |
| `advCoeffx`     | `Real` | advection speed in the x-direction                 | 5e-4     |
| `advCoeffy`     | `Real` | advection speed in the y-direction                 | 5e-4     |
| `diffCoeffx`    | `Real` | diffusion coefficient in the x-direction           | 2e-5     |
| `diffCoeffy`    | `Real` | diffusion coefficient in the y-direction           | 2e-5     |
