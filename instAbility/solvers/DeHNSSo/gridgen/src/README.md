# gridgen/src — layout

Sources are grouped by role in the grid-generation pipeline. `main_gridgen.m`
(one level up) orchestrates them; add the tree to the path with
`addpath(genpath('gridgen/src'))`.

## Pipeline (in call order)

```
load_bf ─► classify_input ─► resample_{cartesian,curvilinear,unstructured} ─► build_stab_grid ─► StabGrid
   (io)          (io)                     (resample)                              (grid)
```

## Folders

| folder | what lives here |
|--------|-----------------|
| `io/` | get the base flow in: `load_bf`, `classify_input` (cartesian / curvilinear / unstructured), `ensure_metadata` (Re/lref/Uref/nu) |
| `resample/` | put the BF on the target grid: `resample_cartesian`, `resample_curvilinear`, `resample_unstructured`, `solve_elliptic_mesh` (curvilinear smoother), `cluster_xi` (streamwise clustering), `report_clustering` (print summary) |
| `grid/` | build the DeHNSSo grid + operators: `build_stab_grid`, `compute_cheb_operators`, `compute_fd_operators`, `compute_bf_derivatives`, `bf_deriv_fields` |
| `forcing/` | wall forcing: `prepare_forcing` (actuator or v-distribution → refinement extent + profile) and `flatten_forcing` (normalise one/many specs); stored as `StabGrid.forcing`, applied by `src/boundary/add_wall_forcing` |
| `math/` | low-level numerics: `deriv1_uneven`, `deriv2_uneven` (non-uniform finite-difference weights), `chebDiff` (Chebyshev) |
| `util/` | small helpers: `get_opt` (param default), `resolve_trim` (fraction → station count) |
| `plot/` | diagnostics: `plot_raw_bf`, `plot_resampled_bf` |

## Note

`math/` (`deriv1_uneven`, `deriv2_uneven`, `chebDiff`) and `grid/build_stab_grid` are
also used by the solver (`src/init/init.m`), which is why solver callers add
`genpath('gridgen/src')` to the path.
