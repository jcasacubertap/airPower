# DeHNSSo — Developer Guide

**DeHNSSo-v2** | Delft Harmonic Navier-Stokes Solver

---

## Contents

- [1. Repository layout](#1-repository-layout)
- [2. Data structures](#2-data-structures)
    - [StabGrid](#stabgrid)
    - [Stab](#stab)
    - [Opt](#opt)
    - [Mode](#mode)
    - [StabRes](#stabres)
    - [Bound](#bound)
- [3. Call flow](#3-call-flow)
- [4. LHS matrix structure](#4-lhs-matrix-structure)
- [5. RHS and inhomogeneous BCs](#5-rhs-and-inhomogeneous-bcs)
- [6. LU strategies](#6-lu-strategies)
- [6a. `lapack_band` MEX](#6a-lapack_band-mex)
- [7. ILST inflow computation](#7-ilst-inflow-computation)
- [8. Adding a new inflow method](#8-adding-a-new-inflow-method)
- [9. DeHNSSo-v2 changes relative to the heritage version](#9-dehnsso-v2-changes-relative-to-the-heritage-version)
- [10. Reference cases](#10-reference-cases)

## 1. Repository layout

```
DeHNSSo/
├── callers/
│   └── benchmark/                  — Per-case stability callers (identical structure)
│       ├── ts_blasius.m            — TS waves in Blasius BL (DeHNSSo-v1 paper Sec. 6.1)
│       ├── ts_compressible.m       — Compressible TS waves, Ma = 0.5 (DeHNSSo-v2 paper Sec. 6.1)
│       ├── cfi_flat_plate.m        — Stationary CFI on flat-plate swept BL (DeHNSSo-v1 paper Sec. 6.2)
│       ├── cfi_hump.m              — Stationary CFI over a hump on swept BL (DeHNSSo-v1 paper Sec. 6.3)
│       ├── cfi_curved_airfoil.m    — Stationary CFI on M3J curved swept-wing BL (DeHNSSo-v2 paper Sec. 6.3)
│       ├── resonator.m             — TS wave over a Helmholtz resonator (DeHNSSo-v2 paper Sec. 6.4)
│       ├── cfi_step.m              — Stationary CFI over a forward-facing step (DeHNSSo-v1 paper Sec. 6.4)
│       ├── wall_jet.m              — Prescribed wall forcing (synthetic jet)
│       └── newton_continuation.m   — Newton pseudo-arclength continuation demo
├── input/                          — DeHNSSo_input_<case>.mat (gridgen → solver handoff, gitignored)
├── output/                         — Saved StabRes results (Opt.save_output, gitignored)
├── literature/                     — Committed reference data (auto-loaded by callers/benchmark/<case>.m via ref_file)
│   ├── DeHNSSo-v1/                 — Three DeHNSSo-v1 cases on v1's own grids and settings, with v1 reference data
│   ├── StabRes_Blasius_NPSE.mat        — NPSE for Blasius TS (Fig. 8)
│   ├── StabRes_SweptWing_flat_NPSE.mat — NPSE for swept-wing flat CFI (Fig. 9)
│   ├── StabRes_SweptWing_Hump_AHLNS.mat — AHLNS for swept-wing hump CFI (Fig. 10)
│   ├── StabRes_SweptWing_step_DNS.mat  — DNS for swept-wing forward-facing-step CFI
│   └── Zanus_PSE_CPG_M05_F40.mat       — PSE for compressible TS, Ma = 0.5
├── src/
│   ├── main.m                      — Main entry point (convergence loop)
│   ├── init/                       — Initialization
│   │   ├── init.m                  — Defaults, mode numbering, grid setup
│   │   ├── inflow.m                — Inflow dispatch + per-mode seeding (local fns ilst_inflow/clst_inflow)
│   │   ├── EVfilter.m              — Physical-mode selection (shared by LST and PSE)
│   │   ├── lst/                    — Local stability theory eigen-solvers
│   │   │   ├── solver_ILST.m       — Incompressible Orr-Sommerfeld/Squire eigenproblem (shared with PSE)
│   │   │   └── solver_CLST.m       — Compressible LST eigenproblem (assembly + polyeig; CLST_Matrix local)
│   │   └── pse/                    — Parabolized stability equations inflow
│   │       ├── solver_IPSE.m       — Incompressible linear PSE (folds in the PSE marcher; seeded by solver_ILST)
│   │       └── solver_CPSE.m       — Compressible linear PSE (Chang & Malik)
│   ├── boundary/                   — Boundary conditions
│   │   ├── setup_boundary.m        — Interpolate Stab.bcw/bct onto the grid; outflow buffer; resolve Opt.hr
│   │   ├── add_wall_forcing.m      — Apply the gridgen-defined wall forcing (sets Stab.bcw)
│   │   ├── apply_resonator_wall.m  — Helmholtz-resonator impedance wall BC + optional mouth u'-slip (Opt.hr)
│   │   └── bc_rows.m               — Row indices of the wall and freestream BC rows in a station block (shared by LHS and RHS)
│   ├── hns/                        — Core HNS equations + solve
│   │   ├── LHS.m                   — Assemble incompressible LHS (per mode)
│   │   ├── LHS_compressible.m      — Assemble compressible (u,v,w,rho,T) LHS
│   │   ├── RHS.m                   — Assemble RHS vector (inflow IC + actuation BCs)
│   │   ├── NLT_HNS.m               — Nonlinear forcing (mode-to-mode interaction)
│   │   └── solve_modes.m           — Mode-by-mode solve, dispatch on Opt.lu_mode
│   ├── linsolve/                   — Direct linear-solver backends (Opt.lu_mode)
│   │   ├── lapack_band/            — LAPACK banded LU MEX (mexbuild_dgbsv.m, dgbsv_*.c)
│   │   └── block_thomas/           — dual-Schur block-Thomas, single-file dispatch (in-core + out-of-core)
│   ├── newton/                     — Newton solver (Opt.solver='newton')
│   │   ├── solve_newton.m          — Pseudo-arclength continuation in A0_fund.
│   │   │                             The only Newton solver: a fixed-amplitude
│   │   │                             solve does not converge at high amplitude.
│   │   ├── newton_context.m        — Shared context: LHS, per-mode LU, BC vectors
│   │   ├── nlt_jacobian.m          — Analytical Jacobian of the NLT forcing (Jh, Ja)
│   │   ├── newton_residual.m       — F(q) in real form [real(q_phys); imag(q_phys)]
│   │   ├── jacobian_apply.m        — Matrix-free Jacobian-vector product
│   │   └── mumps_solve.m           — MUMPS back-end: factor / bsolve / destroy
│   ├── convergence/
│   │   ├── check_convergence.m     — Convergence criterion, amplitude ramping
│   │   └── update_modes.m          — Update u/v/w/p, amplitudes, active mode set
│   ├── io/
│   │   ├── postprocess.m           — Normalize shape functions, compute wavenumbers, StabRes.E
│   │   └── print_config.m          — Print the run configuration
│   ├── plot/
│   │   ├── plot_amplitudes.m       — Final amplitude vs x (with optional reference)
│   │   ├── plot_stability.m        — N-factor, growth rates, wavenumber, phase speed, eigenfunctions
│   │   ├── plot_intermediate.m     — Amplitude during iteration (Figure 2)
│   │   └── plot_convergence.m      — Convergence history (Figure 3, log scale)
│   └── math/
│       ├── chebdif.m              — Chebyshev differentiation matrices (Weideman-Reddy)
│       ├── MappingMalik.m         — Malik wall-normal grid mapping (used by the inflow solvers)
│       └── FD*.m / fornberg_weights.m — Finite-difference operators (uniform + arbitrary spacing)
│   ├── thermo/                     — sutherland.m (viscosity and conductivity laws)
│   └── external/mumps_5.3.5/       — MUMPS sources, MATLAB gateway and install scripts (see docs/dependencies.md)
├── gridgen/                        — Grid generation and base-flow preprocessing
│   ├── main_gridgen.m              — Pipeline entry point (function called by per-case scripts)
│   ├── benchmark/                  — Per-case grid-gen scripts
│   │   ├── one per case, named as in callers/benchmark/: ts_blasius.m, ts_compressible.m, cfi_flat_plate.m,
│   │   │   cfi_hump.m, cfi_step.m, cfi_curved_airfoil.m, resonator.m, wall_jet.m, newton_continuation.m
│   └── src/                        — grouped by pipeline role (see gridgen/src/README.md)
│       ├── io/                     — load_bf, classify_input, ensure_metadata, detect_wall_solid
│       ├── resample/              — resample_{cartesian,curvilinear,unstructured}, solve_elliptic_mesh, cluster_xi, report_clustering
│       ├── grid/                   — build_stab_grid, step_geometry, compute_{cheb,fd}_operators, compute_bf_derivatives, bf_deriv_fields
│       ├── forcing/                — prepare_forcing, flatten_forcing (wall forcing)
│       ├── math/                   — deriv1_uneven, deriv2_uneven, chebDiff
│       ├── util/                   — get_opt, resolve_trim
│       └── plot/                   — plot_raw_bf, plot_resampled_bf
└── baseflow/                       — Boundary-layer solvers
    ├── caller/                     — caller_IBL, caller_CBL (write bf_<case>.mat)
    ├── IBL/                        — Incompressible BL (Blasius, Falkner-Skan-Cooke)
    ├── CBL/                        — Compressible BL (compressibleFSC, continuity_V)
    ├── math/                       — chebdif, MappingMalik, FD1d4o (copies of src/math, kept separate on purpose)
    └── output/benchmark/           — Committed reference base flows: bf_ts_blasius, bf_ts_compressible, bf_cfi_flat_plate
                                      (the M006 resonator flow and the two CSV flows are local only, not shipped)
```

---

## 2. Data structures

### StabGrid

The stability grid, assembled by the gridgen pipeline and stored in the input `.mat` file.

| Field | Size | Description |
|---|---|---|
| `x` | `1 × nx` | Streamwise Reynolds number stations |
| `xun` | `1 × nx` | Physical x-coordinate stations |
| `y` | `ny × 1` | Wall-normal coordinate (index 1 = freestream, ny = wall) |
| `U, V, W` | `ny × nx` | Base-flow velocity components |
| `dxU, dxV, dxW` | `ny × nx` | Streamwise derivatives of base flow |
| `dyU, dyV, dyW` | `ny × nx` | Wall-normal derivatives of base flow |
| `dxxU/V/W/T, ddyU/V/W/T, dxyU, dxyV` | `ny × nx` | Second `∂²/∂x²`, `∂²/∂y²` and mixed `∂²/∂x∂y` base-flow derivatives (compressible) |
| `rho, T, dxrho, dyrho, dxT, dyT` | `ny × nx` | Density/temperature and their first derivatives (compressible) |
| `mu, dmuT, d2muT2` | `ny × nx` | Sutherland viscosity and its `T`-derivatives (compressible) |
| `dzdy, d2zdy2` | `ny × 1` | Malik mapping derivatives |
| `Re` | scalar | Reynolds number at inflow |
| `nx, ny` | scalar | Grid dimensions |
| `H` | scalar | Domain height |

**Grid convention:** index 1 is the freestream (`y = H`), index `ny` is the wall (`y = 0`). `solver_ILST` rebuilds the same Malik grid from `StabGrid.H`/`StabGrid.y_i` and returns the eigenfunction on `StabGrid.y`.

**Base-flow differentiation lives entirely in gridgen.** Every base-flow derivative the LHS assembly and the inflow solvers (`ILST`, `CLST`) need is precomputed by `gridgen/src/grid/compute_bf_derivatives.m` and stored on `StabGrid` — the solver never applies `D1`/`D2` to the base flow. The canonical field list is `gridgen/src/grid/bf_deriv_fields.m`; any of these fields present on the raw base flow are taken as analytic input and skip finite differencing (see User Guide §3, Step 1d). `main.m` copies the required fields into the `BF` struct and errors clearly if a compressible grid predates this (re-run gridgen).

### Stab

User-specified physics parameters.

| Field | Description |
|---|---|
| `M, N` | Spectral truncation (temporal, spanwise) |
| `omega_0, beta_0` | Fundamental frequency and wavenumber |
| `IC` | Inflow method: `'ILST'` |
| `phaseRef` | Phase reference: `'pwall'` or `'umax'` |
| `A0_fund` | Fundamental amplitude (scalar; expanded to full A0 in init.m) |
| `A0` | Initial amplitude vector (1 × nf); built from A0_fund if not set directly |
| `bcx` | x-locations for inhomogeneous BCs |
| `bcw(k,ix,j)` | Wall BC values (3 × nx_bc × nf) |
| `bct(k,ix,j)` | Freestream BC values (3 × nx_bc × nf) |

### Opt

Solver control options.

| Field | Default | Description |
|---|---|---|
| `bc_top` | `{'H_DR','H_DR','H_DR'}` | Freestream BC per component {u,v,w}: `'H_DR'`, `'IH_DR'`, `'H_NM'`, `'free'` |
| `hr` | `[]` | Helmholtz-resonator impedance wall BC over a mouth (struct: `Y` scalar or profile, or `mag`/`phi`; optional `mask`/`x`/`taper`/`uslip`); `[]` = rigid wall |
| `buffer` | `'on'` | Buffer zone type |
| `xb` | 85 | Buffer start (% of domain) |
| `nltbufxb` | = `xb` | NLT-forcing buffer start (% of domain); begins the nonlinear buffer upstream of the linear buffer `xb` |
| `kappa` | 6 | Buffer ramp steepness (larger = steeper) |
| `Conv` | 1e-4 | Convergence criterion |
| `ConvF` | 100 | Relaxation factor during ramping |
| `relax` | off | Under-relaxation factor ω∈(0,1) for the Picard step; stabilizes oscillating cases (omit or ≥1 = off) |
| `anderson` | 0 | Anderson acceleration depth m (0 = off, plain Picard; typical 3–6) |
| `solver` | `'picard'` | Nonlinear solver: `'picard'` or `'newton'` (pseudo-arclength continuation), chosen up front |
| `AFg` | 1.1 | Amplitude factor growth rate |
| `AMAX` | 0.1 | Maximum amplitude factor |
| `TH` | 1e-11 | NLT activation threshold |
| `lu_mode` | `'full'` | LU strategy: `'full'`, `'auto'`, `'seq'`, `'none'`, `'lapack_band'`, `'block_thomas'` (compressible) |
| `lu_max_cache` | 3 | Max cached LU entries for `'auto'` mode |
| `rerun` | `'off'` | Reuse LHS from previous call |
| `save_output` | `[]` | Path to save final `StabRes, StabGrid, Stab` after post-process |
| `save_checkpoint` | `[]` | Path to save pre-postprocess state (`.q` intact) for `IC='LOAD'` |
| `parfor` | auto | Parallel LHS/solve: `'on'` or `'off'` (auto-detected) |
| `linear` | `'off'` | Linear mode: `'on'` (no NLT, 1 iter) or `'off'` |
| `sym_inflow` | off | Inflow modes built as the 2D symmetry image of a canonical mode (`'on'`, `Wbar=0`) |
| `alpha_shift` | `[]` | Compressible LST inflow: pick the eigenvalue nearest this α (continuation seed within a physical band) to de-spike a sweep |
| `drop_alpha2` | false | Compressible LST inflow: `true` zeroes the α² block (CCxx) in the CLST eigenvalue problem to reproduce the old linear GEP; default keeps the full quadratic-in-α operator |
| `plot` | `'on'` | Intermediate plots: `'on'` or `'off'` |
| `verbose` | `'on'` | Console progress output: `'on'` or `'off'` |
| `mmry` | 0 | Internal: 1 = skip LHS assembly (set by rerun) |
| `AF` | scalar | Current amplitude factor (set internally) |
| `dal` | vector | Current relative error per mode (set internally) |
| `LU_cache` | struct | Cached LU factorizations (set internally) |

### Mode

Internal mode indexing (set by `init.m`).

| Field | Description |
|---|---|
| `nf` | Total number of modes |
| `RunJ` | Indices of active (nonzero) modes |
| `Aactive` | 0/1 flag per mode |
| `Mvec, Nvec` | Frequency and wavenumber index per mode |

### StabRes

Solution arrays.

| Field | Size | Description |
|---|---|---|
| `phi(j,:,i)` | `nf × 4*ny × nx` | Concatenated [u v w p] eigenvector |
| `phiIC(j,:)` | `nf × 4*ny` | Inflow IC (reference, frozen) |
| `u,v,w,p(j,y,i)` | `nf × ny × nx` | Shape functions |
| `alpha(j,i)` | `nf × nx` | Streamwise wavenumber |
| `A(j,i)` | `nf × nx` | Mode amplitude |

### Bound

Interpolated BC arrays (set by `setup_boundary.m`).

| Field | Size | Description |
|---|---|---|
| `bcw(k,i,j)` | `3 × nx × nf` | Wall BC on numerical grid |
| `bct(k,i,j)` | `3 × nx × nf` | Freestream BC on numerical grid |
| `*_buf` | — | Buffer-zone arrays |

---

## 3. Call flow

```
callers/benchmark/<case>.m   (ts_blasius / ts_compressible / cfi_flat_plate /
                              cfi_hump / cfi_curved_airfoil / resonator / wall_jet)
└── main(Stab, StabGrid, Opt [, LHS])
    ├── init(StabGrid, Opt, Stab)          % defaults, mode numbering
    ├── inflow(Mode, Stab, StabGrid, StabRes, Opt)
    │   └── ilst_inflow / clst_inflow  (local seeding fns in inflow.m; or solver_IPSE/solver_CPSE)
    │       └── solver_ILST(Re, U, W, ...)  % Orr-Sommerfeld/Squire eigenproblem
    ├── setup_boundary(StabGrid, Opt, Stab)       % interpolate bcw/bct
    ├── print_config(...)
    └── [convergence loop]
        ├── solve_modes(BF, Bound, ...)
        │   ├── LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt)  % LHS assembly
        │   ├── RHS(Bound, Grid, Mode, Opt, StabRes, j)              % RHS + BCs
        │   └── [LU solve with cache]
        ├── AVAL(...)                        % inflow amplitude damping
        ├── update_modes(...)
        ├── check_convergence(...)
        ├── plot_intermediate(...)
        └── plot_convergence(...)
    └── postprocess(...)
```

---

## 4. LHS matrix structure

The LHS is assembled by `LHS.m` for each mode `j`. The matrix encodes the linearised NS equations on the `(ny × nx)` grid, discretized with Chebyshev collocation in `y` and second-order finite differences in `x`. Layout (row block per streamwise station `i`, column block per station):

- Variables per station: `[u, v, w, p]` each of length `ny` → 4*ny rows/cols per station.
- Total size: `4*ny*nx × 4*ny*nx`.

**Key index conventions:**

```matlab
bc_top     = [1, ny+1, 2*ny+1];   % freestream rows (u,v,w at y=H) — per station
bc_top_mfd = [1, 2*ny+1];         % freestream rows for MFD mode (skip v)
bc_wall    = [ny, 2*ny, 3*ny];    % wall rows (u,v,w at y=0) — per station
```

**Freestream BC** (per component via `Opt.bc_top{k}`):
- `'H_DR'` / `'IH_DR'`: identity rows at `bc_top`. RHS injects `bct` value for `IH_DR`.
- `'H_NM'`: rows from `D1f(1,:)` (Chebyshev derivative row at freestream, evaluating `d/dy = 0`).
- `'free'`: no row replacement — the governing equation is kept at `y = H` (as the MFD does for v), used with the parabolization buffer.

**Wall BC:** always Dirichlet (no-slip), except over a resonator mouth when `Opt.hr` is set. `setup_boundary` resolves the settings once into `Bound.hr` (forming `Y` from `mag`/`phi` and the mouth from `StabGrid.forcing` where needed), and `apply_resonator_wall` replaces the wall `v'`-row with the reactive relation `v' − Y·p' = 0` (compressible: on the `rho'`/`T'` columns via the EOS), applied to `MLA1` only so the `(0,0)` MFD keeps no-slip. `Y` may be one value per mouth node for a resolved profile `Y(x)`. Set `hr.uslip` for an optional free-slip `u'` (`du'/dy = 0`) over the mouth.

**`MFD_diff` compression:** the `(0,0)` mean-flow distortion mode (index `jmfd`) does not depend on `omega`/`beta`. Its LHS is stored as `MFD_diff = MLAMFD - MLA1` so that only the difference needs to be recomputed when the base-flow LHS `MLA1` changes.

---

## 5. RHS and inhomogeneous BCs

`RHS.m` builds the right-hand side vector for mode `j`:

1. **Inflow IC:** `R(1:4*ny) = AF * phiIC(j,:)` — the inflow eigenfunction, scaled by the current amplitude factor.
2. **Wall and freestream actuation:** for each station `i` and each velocity component `k=1,2,3`:
   - Freestream: `R(base + (k-1)*ny + 1) = bct(k,i,j)` (first row in the block = y=H row).
   - Wall: `R(base + k*ny) = bcw(k,i,j)` (last row in the block = y=0 row).

Freestream values are only injected for components where `Opt.bc_top{k} = 'IH_DR'`. For `'H_DR'` and `'H_NM'` components, the freestream RHS stays zero. Wall values are always injected (wall BCs are always Dirichlet).

---

## 6. LU strategies

`solve_modes.m` groups the modes by their operator `ML(omega_j, beta_j)`, factors each group once and solves its modes against that factor; `Opt.lu_mode` chooses the factorization and how many factors are kept. Optional iterative refinement (`Opt.iterref = n`) adds `n` passes on any back-end that holds a factor, recovering accuracy on the ill-conditioned operator.

**Sparse LU** (`'full'`, `'auto'`, `'seq'`, `'none'`): `do_lu` runs MATLAB's sparse `lu` (symrcm reordering on the compressible operator, COLAMD otherwise).

- `'full'`: cache every factor in `Opt.LU_cache`, keyed by `(omega_j, beta_j)`. Fastest for repeated solves; memory scales with `nf × (4·ny·nx)² / 8`.
- `'auto'`: LRU cache capped at `Opt.lu_max_cache`. Eviction skips keys still needed this iteration; the cache grows past `lu_max_cache` if every active key is needed (one-time warning).
- `'seq'`: one factor at a time — factor a key, solve every mode that uses it, then discard. No cache, minimal memory.
- `'none'`: no reuse, plain `\` each time. Slowest, minimal memory.

**Banded LU** (`'lapack_band'`): LAPACK banded factor/solve via MEX (`src/linsolve/lapack_band/`; `dgbsv_factor.c` → `zgbtrf`, `dgbsv_solve_cached.c` → `zgbtrs`). Exploits the block-pentadiagonal structure with Apple Accelerate / OpenBLAS dense-banded BLAS3. Factors once per `(omega_j, beta_j)` and caches the factor (LRU, `Opt.lu_max_cache`), so Picard iterations reuse it. ~3× faster than the sparse path on a single solve. See section 6a.

**Block-Thomas** (`'block_thomas'`, compressible only): dual-Schur block-pentadiagonal factorization (Paredes 2024, `src/linsolve/block_thomas/`), partitioned by `Opt.bt_np` (default 2). Factors once per `(omega_j, beta_j)` and caches like the banded path; `Opt.lu_cache='off'` forces a one-shot factor+solve, and `Opt.bt_ooc='on'` streams the factor to disk (out-of-core) for minimal RAM, under `Opt.bt_ooc_dir` (default `tempdir`). `Opt.bt_apply` selects whether each diagonal station block is stored as its pivoted LU factors (default) or as its inverse ('inverse', Paredes 2024); the former is three orders more accurate on near-singular station blocks for ~4% more time; both forms are streamed out of core. Errors on the incompressible operator, whose zero-diagonal continuity rows make the unpivoted block factorization silently wrong.

**MUMPS** (`src/newton/mumps_solve.m`, Newton only): factorises the real 2N x 2N Jacobian through one entry point, `mumps_solve('factor'|'bsolve'|'destroy')`. Not an `Opt.lu_mode` back-end: those factorise one mode at a time, whereas the Newton Jacobian couples all modes at once. It climbs the `Opt.newton_mumps_relax` ladder (ICNTL(14)) before falling back to out-of-core (ICNTL(22)), and deletes its own scratch on destroy, since MUMPS closes those files without unlinking them.

**Rerun:** when `Opt.rerun = 'on'`, the factor cache from the previous run is passed back in via `LHS.LU_cache`. If `(omega_j, beta_j)` is unchanged (e.g. an amplitude sweep) the factor is reused at no cost; if `omega_j` changes (frequency sweep) that mode is refactored. This applies to every caching mode (`'full'`, `'auto'`, `'lapack_band'`, `'block_thomas'`); `'seq'` and `'none'` never cache.

**Rerun with `Opt.hr`:** the resonator rows live in `MLA1`, so a cached factor belongs to one admittance. `main` records what the cached `MLA1` carries (`LHS.hr_sig`, `LHS.hr_rows`) and drops `LHS.LU_cache` whenever it differs. The assembly survives regardless: `solve_modes` recovers `MLAMFD = MLA1 + MFD_diff` (invariant under the wall rows), resets any rows a previous resonator wrote, then applies the current one — so any transition reuses the assembly, including back to a rigid wall. `Opt.hr_in_MLA1` keeps Picard from rewriting the rows every iteration. `newton_context` (Newton) always assembles fresh and so needs no reset.

---

## 6a. `lapack_band` MEX

`Opt.lu_mode = 'lapack_band'` calls a MEX wrapper for LAPACK's banded LU with
partial pivoting. The wrapper auto-detects real vs complex `ML` and dispatches:

- `mxIsComplex(ML) || mxIsComplex(RHS)` → `zgbsv` (complex band LU)
- otherwise → `dgbsv` (real band LU)

The second condition matters in nonlinear runs: the MFD (0,0) mode has a real
`ML` but the NLT-projected RHS may be stored as complex. The MEX promotes
`ML` to complex on the fly when only the RHS is complex.

**Bandwidth detection.** `KL` (subdiagonals) and `KU` (superdiagonals) are
inferred from the sparsity pattern at runtime. For TS at nx=800, ny=81 (4-var):
`KL≈440, KU≈600, LDAB≈1481`, AB storage 6.1 GB complex / 3.1 GB real.

**Why pivoting matters.** The shifted-forward FD boundary stencil at i=2
produces a small diagonal block (`rcond ≈ 1e-20`). Unpivoted block-LU (the
removed `block_penta_solve`) amplified the small pivot through Schur updates
and lost all precision within 4-5 steps. `zgbsv`'s partial pivoting selects
the largest-magnitude candidate within the column band, swaps it to the
diagonal, and proceeds stably. This is why `lapack_band` works on the same
matrices `block_penta_solve` failed on.

**Files:**

- `src/linsolve/lapack_band/dgbsv_solve.c` — MEX source (~150 lines).
- `src/linsolve/lapack_band/mexbuild_dgbsv.m` — one-line build wrapper invoking `mex -R2018a … -lmwlapack`.
- The compiled `.mex*` binary is platform-specific and not version-controlled.

**Build requirements:** MATLAB R2018a+ (for `-R2018a` interleaved-complex API),
a C compiler configured via `mex -setup` (Xcode CLT on macOS, gcc on Linux,
MSVC/MinGW on Windows). MATLAB's LAPACK (`lapack.h` + `libmwlapack`) is used
for portability — on macOS this resolves through Apple Accelerate, on Linux
through OpenBLAS/MKL.

**Factor caching:** each operator is factored once (`dgbsv_factor` -> `*gbtrf`)
and the banded factor is cached in `Opt.LU_cache`, then reused across Picard
iterations (`dgbsv_solve_cached` -> `*gbtrs`), so a nonlinear run pays the
`O(N*KL*KU)` factorization only once per mode. The bands are large, so the cache
is capped by `Opt.lu_max_cache` with LRU eviction (keys still needed by the
current iteration are never evicted). Set `Opt.lu_cache = 'off'` (or
`Opt.lu_max_cache <= 0`) for a one-shot factor+solve that keeps only one band
resident, at the cost of re-factoring on every solve.

---

## 7. ILST inflow computation

`inflow.m` seeds each mode via its local function `ilst_inflow`, which calls
`solver_ILST.m` at the inflow station to find the local stability eigenvalue and
eigenfunction. `solver_ILST` is shared with the incompressible PSE inflow
(`solver_IPSE`, which folds in the PSE marcher); the compressible path is
`clst_inflow` → `solver_CLST.m`. The
eigen-solvers live in `init/lst/`; mode selection, normalization, and storage
into `StabRes` stay in `inflow.m`.

**Key details:**
- `solver_ILST` solves the quadratic (streamwise-viscous `alpha^2` retained) Orr-Sommerfeld/Squire problem; it rebuilds the Malik grid from `StabGrid.H`/`y_i` and returns the eigenfunction on `StabGrid.y`.
- The eigenvalue is filtered by `EVfilter.m` based on the physical growth rate and spatial shape.
- The eigenfunction is normalized so that `max(|u|) = 1` (using spline interpolation for accuracy).
- The phase is then set using `Stab.phaseRef`:
  - `'pwall'`: phase of wall pressure `p(end)` is zeroed.
  - `'umax'`: phase at the `u`-maximum is zeroed.
- Conjugate modes are set by complex conjugate symmetry.

The seed mode decays to zero at the freestream, so the freestream BC type is
immaterial for the inflow eigenproblem (`solver_ILST` clamps it); `Opt.bc_top`
still governs the freestream BC of the HNS march itself.

---

## 8. Adding a new inflow method

To add a method beyond the built-in `'ILST'`/`'CLST'`/`'PSE'`:

1. Add the eigen-solver as a new file in `src/init/lst/` (or `src/init/pse/`) that returns the spectrum and eigenfunctions.
2. Add an `elseif Stab.IC=="..."` branch in `inflow.m` (or a local seeding function like `ilst_inflow`) that calls it, selects the mode via `EVfilter`, and normalizes by `max|u|`.
3. Populate `StabRes.phi(j,:,1)`, `StabRes.u/v/w/p(j,:,1)` at inflow and append active mode indices to `nonzeromodesic`.

---

## 9. DeHNSSo-v2 changes relative to the heritage version

The heritage version is the published DeHNSSo monolithic function (Westerbeek et al., CPC 2024).

The table below is not exhaustive: the compressible operator, the curvilinear formulation, the Newton pseudo-arclength continuation and the direct-solver back-ends are listed in `CHANGELOG.md`.

| Feature | Heritage | DeHNSSo-v2 |
|---|---|---|
| Code structure | Single 900-line function | Modular: ~15 focused files |
| LHS assembly | 6 redundant 14-line BC loops | Vectorized index vectors, no duplication |
| LHS reuse | None | `Opt.rerun`: pass LHS across calls |
| LU cache | Pre-factorized `dML` passed in | `Opt.LU_cache` keyed by `(omega,beta)` |
| Inhomogeneous BCs | Not supported | `Stab.bcw`, `Stab.bct` interpolated and injected in RHS |
| Freestream BC type | Dirichlet only | `Opt.bc_top`: per-component `{'H_DR','IH_DR','H_NM','free'}` |
| ILST phase reference | Fixed (wall pressure) | `Stab.phaseRef`: `'pwall'` or `'umax'` |
| Convergence plot | Linear y-axis (bug) | Log y-axis fixed |
| Caller | Flat script | Sectioned with inline documentation |
| Parallel LHS | None | `Opt.parfor = 'on'`: parfor over stations (LHS) and modes (solve) |
| Intermediate plots | Always on | `Opt.plot = 'on'`/`'off'`; suppressed by default in benchmarks |
| Linear mode | Not supported | `Opt.linear = 'on'`: single iteration, no NLT (linear PSE) |
| Oblique mode | Not supported | `N = 0` with `beta_0 > 0` for single oblique wave |
| Amplitude interface | Manual `ModeToModeNumber` | `Stab.A0_fund` scalar; conjugates built in `init.m` |
| Stability plots | Not available | `plot_stability(StabGrid, StabRes, Opt)`: N-factor, growth rates, eigenfunctions |

---

## 10. Reference cases

Per-case stability runs live in `callers/benchmark/`. The scripts share an identical layout (Setup → Spectral content → Inflow → BCs → Buffer → Solver options → Performance → Nonlinear solver choice → Adjoint → Run → Post-processing) and differ only in case-specific parameter values:

| Case | `Stab.M` | `Stab.N` | `Stab.omega_0` | `Stab.beta_0` | `Stab.phaseRef` | `Opt.buffer` |
|---|---|---|---|---|---|---|
| `blasius` | 1 | 0 | 0.034576 | 0 | `'pwall'` | `'para'` |
| `cfi_flat_plate` | 0 | 1 | 0 | `2π·lref/7.5e-3` | `'umax'` | `'para'` |
| `cfi_hump` | 0 | 1 | 0 | `2π·lref/7.5e-3` | `'umax'` | `'para'` |
| `cfi_curved_airfoil` | 0 | 1 | 0 | 0.40 | `'umax'` | `'on'` |

The corresponding base flows are produced by `gridgen/benchmark/<case>.m` (which writes `input/DeHNSSo_input_<case>.mat`, one file per case so they cannot clobber each other).
