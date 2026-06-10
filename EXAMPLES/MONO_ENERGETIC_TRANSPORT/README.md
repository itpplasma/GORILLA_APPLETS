# Mono-energetic radial transport coefficient

This example computes the mono-energetic radial transport coefficient
D11 with GORILLA and scans it over the normalized collisionality nu*.
It reproduces the transport-coefficient workflow from M. Eder et al.,
*Quasi-geometric integration of guiding-center orbits in piecewise
linear toroidal fields*, Phys. Plasmas 27, 122508 (2020),
<https://doi.org/10.1063/5.0022117>. The reference code for comparison
is NEO-2.

The equilibrium is the bundled AUG g-file `g_file_for_test`
(shot #26884 at 4300 ms), the same axisymmetric tokamak field used by
the other GORILLA examples.

## Run

```bash
cd ../..        # GORILLA_APPLETS root
make            # builds BUILD/gorilla_applets_main.x
cd EXAMPLES/MONO_ENERGETIC_TRANSPORT
./run.sh
```

`run.sh` does both passes in a fresh `run_output/` subdirectory:

1. `i_option=1` precomputes the flux-tube volume `fluxtubevolume.dat`.
   It fixes the Monte Carlo start positions on the chosen flux surface.
2. `i_option=3` traces `n_particles` guiding centers per collisionality
   with the Lorentz operator and fits D11 from the growth of
   `<(Delta s)^2>`. It writes the scan to `nustar_diffcoef_std.dat`.

Set `MODE=single ./run.sh` for one nu* (`i_option=2`) instead of the
scan. Thread count follows `OMP_NUM_THREADS` (default: all cores).

## Output

`run_output/nustar_diffcoef_std.dat` holds one row per collisionality:
nu*, D11, and the standard deviation of the mean. Plot D11 over nu* and
overlay the NEO-2 result for the same configuration.

## Input files

| File | Role |
|------|------|
| `gorilla_applets.inp` | selects the application via `i_option` |
| `mono_energetic_transp_coef.inp` | particle count, energy, collisions, nu* scan |
| `gorilla.inp` | GORILLA core: electron species, flux coordinates, no E-field |
| `tetra_grid.inp` | tetrahedral grid, equilibrium file, grid resolution |
| `seed.inp` | random seed loaded when `seed_option=2`, for reproducible runs |
| `field_divB0.inp`, `preload_for_SYNCH.inp` | g-file field setup (symlinks to core) |
| `MHD_EQUILIBRIA/` | equilibria (symlink to core repo) |

## Resolution and cost

The shipped settings are sized for a quick run: a 30x30x30 grid,
`n_particles=100`, and an 8-point nu* scan. Peak memory is about 18 GB,
dominated by the grid. Raise `n1,n2,n3` in `tetra_grid.inp`,
`n_particles`, and `n_nu_scans` for converged, publication-grade curves;
both memory and runtime grow steeply with grid size.
