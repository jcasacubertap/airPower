# Reference cases

The reference cases below ship with the repository; a base flow marked (not shipped) has to be obtained separately. The first three reproduce published validations in Westerbeek et al., *DeHNSSo: The Delft Harmonic Navier-Stokes Solver for nonlinear stability problems with complex geometric features*, **Comput. Phys. Commun. 302 (2024) 109250** ([DOI](https://doi.org/10.1016/j.cpc.2024.109250)). The others are DeHNSSo-v2 additions, validated against other references.

| Case | Reference | Caller | Base flow |
|---|---|---|---|
| Blasius TS — nonlinear | §6.1, Fig. 8 | `callers/benchmark/ts_blasius.m` | `bf_ts_blasius.mat` |
| Swept-wing CFI flat — nonlinear | §6.2, Fig. 9 | `callers/benchmark/cfi_flat_plate.m` | `bf_cfi_flat_plate.mat` |
| Swept-wing CFI hump — linear | §6.3, Fig. 10 | `callers/benchmark/cfi_hump.m` | `bf_cfi_hump.mat` (not shipped) |
| Compressible TS — linear | Zanus/Bertolotti | `callers/benchmark/ts_compressible.m` | `bf_ts_compressible.mat` |
| Helmholtz resonator — linear | Michelis 2023 | `callers/benchmark/resonator.m` | `bf_compressible_m006_fine.mat` (not shipped) |
| Swept-wing CFI, curved M3J surface | DNS | `callers/benchmark/cfi_curved_airfoil.m` | `bf_cfi_curved_airfoil.csv` (not shipped) |
| Wall jet — prescribed forcing | DeHNSSo-v1 benchmark | `callers/benchmark/wall_jet.m` | `bf_ts_blasius.mat` |
| Newton continuation — demo | — | `callers/benchmark/newton_continuation.m` | `bf_ts_blasius.mat` (400 x 40) |
| Swept-wing CFI step — embedded boundary, linear | DNS | `callers/benchmark/cfi_step.m` | `bf_cfi_step.mat` (not shipped) |

Names are shared by all three stages of a case: `baseflow/output/benchmark/bf_<case>`,
`gridgen/benchmark/<case>.m` writing `input/DeHNSSo_input_<case>.mat`, and
`callers/benchmark/<case>.m`. `ts_` cases are Tollmien–Schlichting, `cfi_` crossflow.

Three of the four cases published with DeHNSSo-v1 (Blasius TS, flat-plate CFI, hump CFI)
are also reproduced on v1's own grids and settings, against the v1 reference data, in
`literature/DeHNSSo-v1/`.

The base flows are committed in `baseflow/output/benchmark/`. Both the gridgen and stability scripts read from there directly, so no manual file copying is required.

> **Shipped configuration.** The unified callers in `callers/benchmark/` ship configured to reproduce the paper figures (nonlinear where the paper is nonlinear), on the medium grid committed in `baseflow/output/benchmark/`. To match a paper figure exactly, only the grid resolution usually needs increasing (regenerate the base flow); each case below notes where.

> **Amplitude convention.** `Stab.A0_fund` is the **physical peak amplitude** `2|û|` of the disturbance — what you set in equals what post-process reports (`StabRes.A`). *(Changed in v2: it was previously the complex amplitude `|û|`, i.e. half this value; see the user guide "Amplitude convention".)*

## Contents

- [Workflow (all cases)](#workflow-all-cases)
- [Case 1 — Blasius TS waves](#case-1--blasius-ts-waves)
- [Case 2 — Swept-wing CFI on a flat plate](#case-2--swept-wing-cfi-on-a-flat-plate)
- [Case 3 — Swept-wing CFI over a smooth hump (linear)](#case-3--swept-wing-cfi-over-a-smooth-hump-linear)
- [Case 4 — Compressible TS waves (Zanus/Bertolotti benchmark)](#case-4--compressible-ts-waves-zanusbertolotti-benchmark)
- [Case 5 — Helmholtz-resonator wall (linear)](#case-5--helmholtz-resonator-wall-linear)
- [Notes for replication](#notes-for-replication)

## Case index

| case | what it is | caller | runtime |
|---|---|---|---|
| [Case 1 — Blasius TS waves](#case-1--blasius-ts-waves) | 2-D nonlinear TS wave packet in a Blasius layer | `ts_blasius` | see Heritage v1 note |
| [Case 2 — Swept-wing CFI on a flat plate](#case-2--swept-wing-cfi-on-a-flat-plate) | nonlinear stationary crossflow, N=5, needs `Opt.anderson = 5` | `cfi_flat_plate` | ~4 days (2544x100, workstation) |
| [Case 3 — Swept-wing CFI over a smooth hump (linear)](#case-3--swept-wing-cfi-over-a-smooth-hump-linear) | linear CFI over a hump, conformal mesh | `cfi_hump` | ~0.6 h (2000x100, desktop) |
| [Case 4 — Compressible TS waves (Zanus/Bertolotti benchmark)](#case-4--compressible-ts-waves-zanusbertolotti-benchmark) | compressible TS, five-variable operator | `ts_compressible` | see Runtime (laptop, single thread) |
| [Case 5 — Helmholtz-resonator wall (linear)](#case-5--helmholtz-resonator-wall-linear) | TS over a resonator mouth, impedance wall BC | `resonator` | ~180 s (4000x40, laptop) |

Each case is documented with the same four subsections: physical parameters, numerical setup, expected output and runtime.

## Workflow (all cases)

```matlab
% 1. Build the stability grid from the base flow
cd gridgen/benchmark/
<case>            % blasius / cfi_flat_plate / cfi_hump
                  % writes ../../input/DeHNSSo_input_<case>.mat

% 2. Run the stability solver
cd ../../callers/benchmark/
<case>            % loads ../../input/DeHNSSo_input_<case>.mat, runs main(), plots
```

`<case>` is the same name in both folders. The stability caller plots amplitude evolution, growth rate, wavenumber, and eigenfunctions.

---

## Case 1 — Blasius TS waves

Tollmien-Schlichting waves in a flat-plate Blasius boundary layer with zero pressure gradient. Reproduces Westerbeek et al. 2024, **Fig. 8** (page 12).

### Physical parameters (paper §6.1)

| Quantity | Value |
|---|---|
| Reference velocity `U_ref` | 10 m/s (constant external velocity) |
| Reference length `δ̄₀` | 6.075×10⁻⁴ m (Blasius length at inflow) |
| Kinematic viscosity `ν` | 1.5188×10⁻⁵ m²/s |
| Reynolds number `Re = δ̄₀·U_ref/ν` | 400 |
| Fundamental frequency `f` | 90.6 Hz |
| Non-dim frequency `ω₁ = 2π·f·δ̄₀/U_ref` | 0.034576 |
| Spanwise wavenumber `β₀` | 0 (2D) |

### Numerical setup

| Knob | Default (shipped) | Paper (Fig. 8) |
|---|---|---|
| Spectral truncation `(M, N)` | `(5, 0)` — 5 harmonics + MFD | same |
| `Opt.linear` | `'off'` (nonlinear) | same |
| `Stab.A0_fund` | `0.0025·√2` (physical peak) | same |
| Streamwise grid `nx` | 800 (medium, from `bf_ts_blasius.mat`) | 2674 (refined) |
| Wall-normal grid `ny` | 40 (medium) | 100 (refined) |
| Buffer | `Opt.buffer = 'on'`, starts at 70% | same |

### To reproduce Fig. 8 exactly (refined grid)

The shipped caller is already the M=5 nonlinear configuration. For the refined paper grid, regenerate `bf_ts_blasius.mat` at `nx=2674`, `ny=100` via `baseflow/caller/caller_IBL.m`.

### Expected output (paper-faithful run)

- Amplitude evolution of modes (1,0) through (5,0) plus the (0,0) mean-flow distortion, matching the NPSE reference within line thickness.
- Visible kink in the MFD curve near `x ≈ 1700` where its peak amplitude switches between two local maxima — physical, present in NPSE too.

### Runtime

- Heritage v1 (Westerbeek 2024, Table A.1)
  - Refined grid (2672×100): ~12 h on a desktop workstation (Intel Xeon, 8 cores, 128 GB)
  - Medium grid (800×40): ~0.45 h on the same hardware
- v2 single-thread on a 32 GB Apple Silicon laptop
  - Linear case (M=1, 800×40): **25.2 s** (heritage 194.7 s → ~7.7× faster)
  - Nonlinear case (M=5, 800×40): **88.5 s** (heritage 785.8 s → ~8.8× faster)
  - LHS stored memory: 115.4 MB (heritage 224.4 MB → 49% less)

---

## Case 2 — Swept-wing CFI on a flat plate

Stationary crossflow instabilities (CFI) in a swept-wing boundary layer with favorable pressure gradient. Reproduces Westerbeek et al. 2024, **Fig. 9** (page 13). Equivalent to *Case A* in Casacuberta et al. (J. Fluid Mech. 2022).

### Physical parameters (paper §6.2)

| Quantity | Value |
|---|---|
| Reference velocity `U_ref` | 15.1 m/s (inflow streamwise component) |
| Reference length `δ̄₀` | 2.14×10⁻⁴ m |
| Kinematic viscosity `ν` | 1.47×10⁻⁵ m²/s |
| Reynolds number `Re` | 220 |
| Spanwise external velocity `W_e` | -1.24 (constant) |
| Streamwise external velocity `U_e(x)` | polynomial fit, eq. (26) of paper |
| Fundamental frequency `ω₀` | 0 (stationary) |
| Spanwise wavenumber `β₀` | `2π·lref / 7.5e-3` (λ = 7.5 mm) ≈ 0.18 |

### Numerical setup

| Knob | Default (shipped) | Paper (Fig. 9) |
|---|---|---|
| Spectral truncation `(M, N)` | `(0, 5)` — 5 harmonics + MFD | same |
| `Opt.linear` | `'off'` (nonlinear) | same |
| `Stab.A0_fund` | `3.5×10⁻³` (physical peak) | same |
| Streamwise grid `nx` | 1200 (after gridgen Cartesian resample) | 1272 |
| Wall-normal grid `ny` | 60 (Chebyshev) | 100 |
| Buffer | `Opt.buffer = 'on'`, starts at 85%; `Opt.anderson = 5` (required for the cascade) | same |
| Amplitude ramp | cold-starts and ramps by `γ = 1.10` per iteration (nonlinear) | same |

### To reproduce Fig. 9 exactly (refined grid)

The shipped caller is already the N=5 nonlinear configuration with `Opt.anderson = 5` — **required**, since without it plain Picard diverges during the mode cascade near saturation (`max_A` runs away). For the paper's wall-normal resolution, increase `params.n_eta_new = 100` in `gridgen/benchmark/cfi_flat_plate.m` and re-run gridgen.

### Expected output (paper-faithful run)

- Amplitude evolution of modes (0,1) through (0,5) plus MFD, matching DNS (INCA) and NPSE references within line thickness across `x ∈ [220, 1400]`.
- Fundamental mode (0,1) saturates at `u'_max ≈ 0.3` (CFI saturation level); secondary lobe in eigenfunction at `x ≈ 1359` (Fig. 9f) is captured.

### Runtime

- Refined grid (2544×100): ~4 days on workstation
- Medium grid (1200×40): ~1.3 h
- v2 single-thread on a 32 GB Apple Silicon laptop: **303 s, 65 Picard iterations, peak 15.1 GB** — ~15× faster than heritage v1 (Westerbeek 2024 Table A.1: 4644 s on 8-core / 32 GB)

---

## Case 3 — Swept-wing CFI over a smooth hump (linear)

Linear development of a stationary CFI mode interacting with a smooth Gaussian hump on the swept-wing surface. Reproduces Westerbeek et al. 2024, **Fig. 10** (page 14).

### Physical parameters (paper §6.3)

Same external flow as Case 2 (`U_e(x)`, `W_e = -1.24`, `Re = 220`, `δ̄₀ = 2.14×10⁻⁴ m`, `U_ref = 15.1 m/s`).

| Hump geometry | Value |
|---|---|
| Hump height `h` | `2.5·δ̄₀` |
| Hump width `b` | `39.9·δ̄₀` |
| Hump centre `x_m` | 859 |
| Wall profile `y_wall(x)` | `h·exp(−((x−x_m)/b)²)` |

The hump is shallow enough to avoid flow separation; the base flow is taken from a steady incompressible NS solution (e.g. COMSOL with second-order velocity elements).

### Numerical setup

| Knob | Default (shipped, matches paper) |
|---|---|
| Spectral truncation `(M, N)` | `(0, 1)` — fundamental + conjugate + MFD (3 modes; MFD inactive in linear mode) |
| `Opt.linear` | `'on'` (linear analysis only) |
| `Stab.A0_fund` | `3.5×10⁻³` (linear regime — value sets only the absolute scale) |
| Streamwise grid `nx` | 1200 (curvilinear, Gaussian-clustered around `x_m = 859`) |
| Wall-normal grid `ny` | 80 (FD, 4th-order, sparse) |
| Buffer | `Opt.buffer = 'para'`, starts at 85% |
| Mesh | grid is **conformal** to the hump (`η`-axis follows the wall) — handled automatically by `gridgen/benchmark/cfi_hump.m` |

The shipped defaults match the paper configuration; no edits needed to reproduce Fig. 10.

### Expected output

- Linear amplitude development (Fig. 10a) showing a slight stabilization downstream of the hump.
- Growth rate `α_i(x)` (Fig. 10c) with the characteristic dip-peak-dip signature around `x ∈ [820, 1000]`, matching the AHLNS reference.

### Runtime

- 2000×100 grid: ~0.6 h on a desktop (Intel Xeon, 4 cores, 32 GB)

---

## Case 4 — Compressible TS waves (Zanus/Bertolotti benchmark)

Tollmien-Schlichting waves in a compressible (subsonic) flat-plate boundary layer at Ma = 0.5. Reproduces the Zanus/Bertolotti (2017) PSE CPG benchmark.

### Physical parameters

| Quantity | Value |
|---|---|
| Mach number `Ma` | 0.50 |
| Prandtl number `Pr` | 0.72 |
| Heat-capacity ratio `γ` | 1.4 |
| Sutherland constant `S* = S/T_∞` | 0.536 (air, `T_∞ = 206 K`) |
| Wall thermal BC | adiabatic |
| Reduced frequency `F = 2π f ν / U_∞²` | 40×10⁻⁶ |
| Inflow Reynolds `R₀ = √Re_x₀` | 500 |
| Outflow Reynolds `R₁` | 1700 |

### Numerical setup

| Knob | Default (shipped) | Paper (M = 0.5) |
|---|---|---|
| Spectral truncation `(M, N)` | `(1, 0)` — linear single TS wave | same |
| `Opt.compressible` | `'on'` | required |
| `Opt.linear` | `'on'` | same |
| `Stab.IC` | `'CLST'` (required for compressible) | same |
| `Opt.bc_top` | `{'H_DR','H_DR','H_DR','H_DR'}` (4-channel u,v,w,T) | same |
| `Opt.bc_wall_thermal` | `'isothermal'` (also OK for adiabatic BF) | same |
| `Opt.buffer` | `'para'` (required at Ma > 0) | same |
| Streamwise grid `nx` | 800 (from compressible FSC) | 800 |
| Wall-normal grid `ny` | 81 (Malik–Chebyshev) | 81 |
| `Opt.lu_mode` | `'lapack_band'` (recommended) or `'seq'` | — |

### Pipeline

```matlab
cd baseflow/caller/;       caller_CBL                          % bf_ts_compressible.mat (~10 s)
cd ../../gridgen/benchmark/; ts_compressible              % DeHNSSo_input_ts_compressible.mat (~10 s)
cd ../../callers/benchmark/; ts_compressible              % HNS + Zanus/Bertolotti overlay
```

If `Opt.lu_mode = 'lapack_band'`: build the MEX once per machine with `run('src/linsolve/lapack_band/mexbuild_dgbsv.m')`. Requires MATLAB ≥ R2018a and a C compiler configured for `mex`. Fall back to `Opt.lu_mode = 'seq'` if no compiler is available (≈3× slower).

### Expected output

- Spatial growth rate `−α_i(R)` matches the Zanus/Bertolotti 2017 PSE reference (markers) to plot resolution across the full domain (R ∈ [500, 1700]), including the unstable lobe peaked at R ≈ 1050 and the return to stability at R ≳ 1400.
- Streamwise wavenumber `α_r(R)` overlays the parallel CLST result.
- `StabRes.E` (Chu/Mack disturbance energy) is written by `postprocess`; `Opt.plot_metric = 'energy'` plots its square root.

### Runtime

v2 single-thread on a 32 GB Apple Silicon laptop:

- `lu_mode = 'lapack_band'`: ~91 s total (13 s LHS + 75 s LU+solve), ~20 GB peak (banded AB)
- `lu_mode = 'seq'`: ~226 s total (13 s LHS + 211 s LU+solve), ~8 GB peak (sparse LU)

---

## Case 5 — Helmholtz-resonator wall (linear)

A Tollmien-Schlichting wave passing over a Helmholtz-resonator neck embedded in the wall, imposed as a reactive impedance boundary condition. Validated against Michelis, de Koning & Kotsonis, *On the interaction of Tollmien-Schlichting waves with a wall-embedded Helmholtz resonator*, **Phys. Fluids 35, 034104 (2023)**.

### Physical parameters

| Quantity | Value |
|---|---|
| Freestream velocity `U_0` | 20.9 m/s |
| Reference length `l_ref` | `2.695×10⁻⁴` m (read from the base flow) |
| Inflow Reynolds number `Re_0` | 371 |
| Mach number | 0.061 (effectively incompressible) |
| TS frequency `f_TS` | 295 Hz, equal to the resonance `f_r` (on resonance) |
| Neck width | 1 mm, centred at `x_0 = 0.55` m |
| Admittance | `\|Y\| = 12.7688`, `arg Y = π`, entered as a constant: the lumped-model value at resonance (the model itself is not part of the tree) |

On resonance the reactances cancel, so the admittance is real and negative: the phase is exactly `π`, which is why a resonant passive resonator can only amplify.

### Numerical setup

| Knob | Default (shipped) |
|---|---|
| Spectral truncation `(M, N)` | `(1, 0)` — single TS mode |
| `Opt.linear` | `'on'` |
| `Opt.compressible` | `'on'` |
| `Stab.IC` | `'CLST'` (compressible local stability, slow mode) |
| `Stab.A0_fund` | `1×10⁻⁴` (arbitrary in linear mode) |
| Streamwise grid `nx` | 4000, Gaussian-clustered at the neck (`spread = 30`, `refine = 5`) |
| Wall-normal grid `ny` | 40 (FD, 4th-order) |
| Buffer | `Opt.buffer = 'on'`, starts at 90% |

The neck geometry lives in `gridgen/benchmark/resonator.m` as `params.forcing`, whose `w` is the **physical mouth half-width**; `params.spread` sets the refinement extent separately. That clustering puts 21 stations across the 1 mm mouth, at the finest spacing on the grid. The caller reads the centre and width back from `StabGrid.forcing`, so the two never drift apart.

### Expected output

```
Resonator wall BC active: |Y| = 12.7688, arg Y = +180.0 deg, mouth = 21 stations.
dN over the flat-plate value: 0.354      (Michelis ~ +0.4, amplifying)
```

Measured 2026-09-11 with the neck centred on x_r = 0.55 m (`Xc = 0.55/lref`; the
earlier `0.55/lref + s_nd/2` put the mouth 0.5 mm downstream of the COMSOL reference
and gave 0.353).

Two solves are run, a rigid-wall baseline and the resonator, and the N-factor increment `ΔN` is their difference just upstream of the buffer. The resonator **amplifies**, consistent with Michelis 2023.

### Runtime

- 4000×40 grid, both solves: ~180 s on a laptop (Apple M-series, 32 GB)

---

## Notes for replication

- **Reference benchmark data.** NPSE (Cases 1, 2) and AHLNS (Case 3) reference curves are shipped in `literature/` and auto-loaded by the unified callers via `ref_file`. To override or skip the overlay, edit the `ref_file = ...` line near the top of the post-processing block in `callers/benchmark/<case>.m`.
- **Buffer treatment.** The hump (Case 3) and compressible (Case 4) cases use parabolisation (`Opt.buffer = 'para'`); the Blasius, swept-wing-flat and resonator (Case 5) cases use amplitude damping (`'on'`). The buffer starts at 70–90% of the streamwise domain depending on the case.
- **Convergence criterion.** `Opt.Conv = 1e-4` (default); for the nonlinear CFI case the absolute floor `1e-14` in `check_convergence.m` is what allows the harmonics to fully converge despite their tiny amplitudes.
- **Memory.** Cases 1, 3 and 5 fit comfortably on a laptop. Case 2 (nonlinear CFI, 6 modes) needs roughly 15 GB of RAM at peak.
