# DeHNSSo — Developer Guide

**Version 1.2** | Delft Harmonic Navier-Stokes Solver

---

## 1. Repository layout

```
DeHNSSo/
├── callers/
│   └── benchmark/                  — Per-case stability callers (identical structure)
│       ├── blasius.m               — TS waves in Blasius BL
│       ├── sweptwing_flat.m        — Stationary CFI on flat-plate swept BL
│       ├── sweptwing_hump.m        — Stationary CFI over a hump on swept BL
│       └── m3j.m                   — Stationary CFI on M3J curved swept-wing BL
├── input/                          — DeHNSSo_input.mat (gridgen → solver handoff, gitignored)
├── literature/                     — Committed reference data (auto-loaded by callers/benchmark/<case>.m via ref_file)
│   ├── StabRes_Blasius_NPSE.mat        — NPSE for Blasius TS (Fig. 8)
│   ├── StabRes_Blasius_LPSE.mat        — LPSE for Blasius TS (linear PSE)
│   ├── StabRes_SweptWing_flat_NPSE.mat — NPSE for swept-wing flat CFI (Fig. 9)
│   └── StabRes_SweptWing_Hump_AHLNS.mat — AHLNS for swept-wing hump CFI (Fig. 10)
├── src/
│   ├── main.m                      — Main entry point (convergence loop)
│   ├── init/                       — Initialisation
│   │   ├── init.m                  — Defaults, mode numbering, grid setup
│   │   ├── inflow.m                — Inflow BC dispatch (ILST or LOAD)
│   │   ├── ILST.m                  — Calls ILST_solver, selects/normalises eigenfunction
│   │   ├── ILST_solver.m           — Orr-Sommerfeld/Squire generalised eigenvalue problem
│   │   └── boundary.m              — Interpolates Stab.bcw/bct onto numerical grid
│   ├── hns/                        — Core HNS equations
│   │   ├── LHS.m                   — Assemble global LHS sparse matrix (per mode)
│   │   ├── RHS.m                   — Assemble RHS vector (inflow IC + actuation BCs)
│   │   ├── solve_modes.m           — Mode-by-mode solve with LU cache
│   │   ├── NLT.m                   — Nonlinear forcing (mode-to-mode interaction)
│   │   └── AVAL.m                  — Inflow amplitude damping
│   ├── convergence/
│   │   ├── check_convergence.m     — Convergence criterion, amplitude ramping
│   │   └── update_modes.m          — Update u/v/w/p, amplitudes, active mode set
│   ├── postprocess/
│   │   ├── postprocess.m           — Normalise shape functions, compute wavenumbers
│   │   └── ModeToModeNumber.m
│   ├── plot/
│   │   ├── plot_amplitudes.m       — Final amplitude vs x (with optional reference)
│   │   ├── plot_stability.m        — N-factor, growth rates, wavenumber, phase speed, eigenfunctions
│   │   ├── plot_intermediate.m     — Amplitude during iteration (Figure 2)
│   │   └── plot_convergence.m      — Convergence history (Figure 3, log scale)
│   └── math/
│       ├── chebDiff.m              — Chebyshev differentiation matrix
│       └── FD*.m                   — Finite-difference operators (2nd/4th order)
├── gridgen/                        — Grid generation and base-flow preprocessing
│   ├── main_gridgen.m              — Pipeline entry point (function called by per-case scripts)
│   ├── benchmark/                  — Per-case grid-gen scripts (blasius, sweptwing_flat, ...)
│   └── src/                        — Mesh generation, metric tensors, scattered/structured resampling
└── baseflow/                       — Boundary-layer solvers
    ├── IBL/                        — Incompressible BL (Blasius, Falkner-Skan-Cooke)
    └── output/benchmark/           — Committed reference base flows (bf_*.mat, bf_m3j.csv)
```

---

## 2. Data structures

### StabGrid

The stability grid, assembled by the Interpol pipeline and stored in the input `.mat` file.

| Field | Size | Description |
|---|---|---|
| `x` | `1 × nx` | Streamwise Reynolds number stations |
| `xun` | `1 × nx` | Physical x-coordinate stations |
| `y` | `ny × 1` | Wall-normal coordinate (index 1 = freestream, ny = wall) |
| `U, V, W` | `ny × nx` | Base-flow velocity components |
| `dxU, dxV, dxW` | `ny × nx` | Streamwise derivatives of base flow |
| `dyU, dyV, dyW` | `ny × nx` | Wall-normal derivatives of base flow |
| `dzdy, d2zdy2` | `ny × 1` | Malik mapping derivatives |
| `Re` | scalar | Reynolds number at inflow |
| `nx, ny` | scalar | Grid dimensions |
| `H` | scalar | Domain height |

**Grid convention:** index 1 is the freestream (`y = H`), index `ny` is the wall (`y = 0`). ILST_solver internally flips the arrays (it expects index 1 = wall).

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
| `bc_top` | `{'H_DR','H_DR','H_DR'}` | Freestream BC per component {u,v,w}: `'H_DR'`, `'IH_DR'`, `'H_NM'` |
| `buffer` | `'on'` | Buffer zone type |
| `xb` | 85 | Buffer start (% of domain) |
| `Conv` | 1e-4 | Convergence criterion |
| `ConvF` | 100 | Relaxation factor during ramping |
| `AFg` | 1.1 | Amplitude factor growth rate |
| `AMAX` | 0.1 | Maximum amplitude factor |
| `TH` | 1e-11 | NLT activation threshold |
| `lu_mode` | `'full'` | LU cache strategy: `'full'`, `'auto'`, `'none'` |
| `lu_max_cache` | 3 | Max cached LU entries for `'auto'` mode |
| `rerun` | `'off'` | Reuse LHS from previous call |
| `parfor` | auto | Parallel LHS/solve: `'on'` or `'off'` (auto-detected) |
| `gpu` | auto | GPU solve via gpuArray: `'on'` or `'off'` (auto-detected) |
| `linear` | `'off'` | Linear mode: `'on'` (no NLT, 1 iter) or `'off'` |
| `plot` | `'on'` | Intermediate plots: `'on'` or `'off'` |
| `mmry` | 0 | Internal: 1 = skip LHS assembly (set by rerun) |
| `AF` | scalar | Current amplitude factor (set internally) |
| `dal` | vector | Current relative error per mode (set internally) |
| `LU_cache` | struct | Cached LU factorisations (set internally) |

### Mode

Internal mode bookkeeping (set by `init.m`).

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

Interpolated BC arrays (set by `boundary.m`).

| Field | Size | Description |
|---|---|---|
| `bcw(k,i,j)` | `3 × nx × nf` | Wall BC on numerical grid |
| `bct(k,i,j)` | `3 × nx × nf` | Freestream BC on numerical grid |
| `*_buf` | — | Buffer-zone arrays |

---

## 3. Call flow

```
callers/benchmark/<case>.m   (blasius / sweptwing_flat / sweptwing_hump / m3j)
└── main(Stab, StabGrid, Opt [, LHS])
    ├── init(StabGrid, Opt, Stab)          % defaults, mode numbering
    ├── inflow(Mode, Stab, StabGrid, StabRes, Opt)
    │   └── ILST(j, 1, StabGrid, StabRes, bc_top)
    │       └── ILST_solver(Re, U, W, ...)  % eigenvalue problem
    ├── boundary(StabGrid, Opt, Stab)       % interpolate bcw/bct
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

The LHS is assembled by `LHS.m` for each mode `j`. The matrix encodes the linearised NS equations on the `(ny × nx)` grid, discretised with Chebyshev collocation in `y` and second-order finite differences in `x`. Layout (row block per streamwise station `i`, column block per station):

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

**Wall BC:** always Dirichlet (no-slip).

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

## 6. LU cache

`solve_modes.m` caches LU factorisations in `Opt.LU_cache`, keyed by `(omega_j, beta_j)`. On the first solve for a mode, the LU is computed and stored. On subsequent iterations, it is reused if the key matches.

**Strategies (`Opt.lu_mode`):**
- `'full'`: unlimited cache — fastest for repeated solves, memory scales with `nf × (4*ny*nx)^2 / 8`.
- `'auto'`: capped at `Opt.lu_max_cache` entries, oldest evicted (LRU). Use when memory is limited.
- `'none'`: no cache — MATLAB's `\` is called every time. Slowest but minimal memory.

**GPU mode (`Opt.gpu = 'on'`):** the assembled `ML` matrix is stored as `gpuArray(ML)` in the cache instead of a `decomposition` object. MATLAB does not support sparse complex LU decomposition caching on GPU; the matrix is refactorised on every iteration. Only beneficial when ny is large enough that O(ny²) GPU compute exceeds host-GPU transfer cost (roughly ny ≥ 80 for NVIDIA A100).

**Rerun:** when `Opt.rerun = 'on'`, the LU cache from the previous run is passed back in via `LHS.LU_cache`. If `(omega_j, beta_j)` is unchanged (e.g., amplitude sweep), the factorisation is reused with no cost. If `omega_j` changes (frequency sweep), a new factorisation is computed for that mode.

---

## 7. ILST inflow computation

`ILST.m` calls `ILST_solver.m` at the inflow station to find the local stability eigenvalue and eigenfunction.

**Key details:**
- Arrays are flipped before passing to `ILST_solver` (internal convention: index 1 = wall).
- The eigenvalue is filtered by `EVfilter.m` based on the physical growth rate and spatial shape.
- The eigenfunction is normalised so that `max(|u|) = 1` (using spline interpolation for accuracy).
- The phase is then set using `Stab.phaseRef`:
  - `'pwall'`: phase of wall pressure `p(end)` is zeroed.
  - `'umax'`: phase at the `u`-maximum is zeroed.
- Conjugate modes are set by complex conjugate symmetry.

**Neumann freestream BC in ILST_solver:**
When any component uses `'H_NM'`, the ILST solver receives `bc_top = 'H_NM'` and the farfield rows `n4-2:n4` use rows from the Chebyshev differentiation matrix `Total_Dn` rather than identity rows, imposing `d/dy = 0` at `y = H`.

---

## 8. Adding a new inflow method

To add a method beyond `'ILST'`:

1. Add a new `elseif` branch in `inflow.m` (around line 76).
2. Populate `StabRes.phi(j,:,1)`, `StabRes.u/v/w/p(j,:,1)` at inflow.
3. Append active mode indices to `nonzeromodesic`.

---

## 9. v1.2 changes relative to heritage version

The heritage version is the published DeHNSSo monolithic function (Westerbeek et al., CPC 2024).

| Feature | Heritage | v1.2 |
|---|---|---|
| Code structure | Single 900-line function | Modular: ~15 focused files |
| LHS assembly | 6 redundant 14-line BC loops | Vectorised index vectors, no duplication |
| LHS reuse | None | `Opt.rerun`: pass LHS across calls |
| LU cache | Pre-factorised `dML` passed in | `Opt.LU_cache` keyed by `(omega,beta)` |
| Inhomogeneous BCs | Not supported | `Stab.bcw`, `Stab.bct` interpolated and injected in RHS |
| Freestream BC type | Dirichlet only | `Opt.bc_top`: per-component `{'H_DR','IH_DR','H_NM'}` |
| ILST phase reference | Fixed (wall pressure) | `Stab.phaseRef`: `'pwall'` or `'umax'` |
| Convergence plot | Linear y-axis (bug) | Log y-axis fixed |
| Caller | Flat script | Sectioned with inline documentation |
| Parallel LHS | None | `Opt.parfor = 'on'`: parfor over stations (LHS) and modes (solve) |
| GPU solve | None | `Opt.gpu = 'on'`: gpuArray `\` for large ny; auto-detected |
| Intermediate plots | Always on | `Opt.plot = 'on'`/`'off'`; suppressed by default in benchmarks |
| Linear mode | Not supported | `Opt.linear = 'on'`: single iteration, no NLT (linear PSE) |
| Oblique mode | Not supported | `N = 0` with `beta_0 > 0` for single oblique wave |
| Amplitude interface | Manual `ModeToModeNumber` | `Stab.A0_fund` scalar; conjugates built in `init.m` |
| Stability plots | Not available | `plot_stability(StabGrid, StabRes, Opt)`: N-factor, growth rates, eigenfunctions |

---

## 10. Reference cases

Per-case stability runs live in `callers/benchmark/`. All four scripts share an identical layout (Setup → Spectral content → Inflow → BCs → Buffer → Solver options → Performance → Nonlinear solver choice → Adjoint → Run → Post-processing) and differ only in case-specific parameter values:

| Case | `Stab.M` | `Stab.N` | `Stab.omega_0` | `Stab.beta_0` | `Stab.phaseRef` | `Opt.buffer` |
|---|---|---|---|---|---|---|
| `blasius` | 1 | 0 | 0.034576 | 0 | `'pwall'` | `'para'` |
| `sweptwing_flat` | 0 | 1 | 0 | `2π·lref/7.5e-3` | `'umax'` | `'para'` |
| `sweptwing_hump` | 0 | 1 | 0 | `2π·lref/7.5e-3` | `'umax'` | `'para'` |
| `m3j` | 0 | 1 | 0 | 0.40 | `'umax'` | `'on'` |

The corresponding base flows are produced by `gridgen/benchmark/<case>.m` (which writes `input/DeHNSSo_input.mat`).
