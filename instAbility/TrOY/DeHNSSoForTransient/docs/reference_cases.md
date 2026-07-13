# Reference cases

Three reference cases ship with the repository. Each reproduces a published validation in Westerbeek et al., *DeHNSSo: The Delft Harmonic Navier-Stokes Solver for nonlinear stability problems with complex geometric features*, **Comput. Phys. Commun. 302 (2024) 109250** ([DOI](https://doi.org/10.1016/j.cpc.2024.109250)).

| Case | Paper section | Paper figure | Caller | Base flow |
|---|---|---|---|---|
| Blasius TS — nonlinear | §6.1 | Fig. 8 | `callers/benchmark/blasius.m` | `bf_blasius.mat` |
| Swept-wing CFI flat — nonlinear | §6.2 | Fig. 9 | `callers/benchmark/sweptwing_flat.m` | `bf_sweptwing_flat.mat` |
| Swept-wing CFI hump — linear | §6.3 | Fig. 10 | `callers/benchmark/sweptwing_hump.m` | `bf_sweptwing_hump.mat` |

The base flows are committed in `baseflow/output/benchmark/`. Both the gridgen and stability scripts read from there directly, so no manual file copying is required.

> **Defaults vs paper-faithful.** The unified callers in `callers/benchmark/` ship with **fast smoke** defaults (linear, single fundamental mode, medium grid). To reproduce the Westerbeek 2024 paper figures exactly, edit the spectral content (`Stab.M`, `Stab.N`), set `Opt.linear = 'off'`, and where noted regenerate the base flow at finer resolution. Each case below has a "To reproduce Fig. X" callout listing the specific changes.

> **Amplitude convention.** `Stab.A0_fund` is the amplitude of the **fundamental mode only**; `init.m` mirrors it onto the conjugate. The peak `u'` in the paper therefore corresponds to `2 · Stab.A0_fund` for a single-fundamental case.

## Workflow (all three cases)

```matlab
% 1. Build the stability grid from the base flow
cd gridgen/benchmark/
<case>            % blasius / sweptwing_flat / sweptwing_hump
                  % writes ../../input/DeHNSSo_input.mat

% 2. Run the stability solver
cd ../../callers/benchmark/
<case>            % loads ../../input/DeHNSSo_input.mat, runs main(), plots
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
| Spectral truncation `(M, N)` | `(1, 0)` — linear single TS wave | `(5, 0)` — 5 harmonics + MFD |
| `Opt.linear` | `'on'` | `'off'` |
| `Stab.A0_fund` | `0.00125·√2` (paper A0 = `2·A0_fund` = `0.0025·√2`) | `0.00125·√2` |
| Streamwise grid `nx` | 800 (medium, from `bf_blasius.mat`) | 2674 (refined) or 800 (medium) |
| Wall-normal grid `ny` | 40 (medium) | 100 (refined) or 40 (medium) |
| Buffer | `Opt.buffer = 'para'`, starts at 85% | same |

### To reproduce Fig. 8 (M=5 nonlinear, medium grid)

In `callers/benchmark/blasius.m`: set `Stab.M = 5` and `Opt.linear = 'off'`. For the refined paper grid, regenerate `bf_blasius.mat` at `nx=2674`, `ny=100` via `baseflow/IBL/caller.m`.

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
| Spectral truncation `(M, N)` | `(0, 1)` — fundamental + conjugate + MFD | `(0, 5)` — 5 harmonics + MFD |
| `Opt.linear` | `'on'` | `'off'` |
| `Stab.A0_fund` | `1.75×10⁻³` (= `3.5e-3 / 2`) | `1.75×10⁻³` |
| Streamwise grid `nx` | 1200 (after gridgen Cartesian resample) | 1272 |
| Wall-normal grid `ny` | 60 (Chebyshev) | 100 |
| Buffer | `Opt.buffer = 'para'`, starts at 85% | same |
| Amplitude ramp | `Opt.linear = 'on'` skips it; for nonlinear runs, cold-start at `A0 ≈ 2.6×10⁻⁴` and ramp by `γ = 1.10` per iteration |

### To reproduce Fig. 9 (N=5 nonlinear)

In `callers/benchmark/sweptwing_flat.m`: set `Stab.N = 5` and `Opt.linear = 'off'`. For the paper's wall-normal resolution, increase `params.n_eta_new = 100` in `gridgen/benchmark/sweptwing_flat.m` and re-run gridgen.

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
| `Stab.A0_fund` | `1.75×10⁻³` (linear regime — value sets only the absolute scale) |
| Streamwise grid `nx` | 1200 (curvilinear, Gaussian-clustered around `x_m = 859`) |
| Wall-normal grid `ny` | 80 (FD, 4th-order, sparse) |
| Buffer | `Opt.buffer = 'para'`, starts at 85% |
| Mesh | grid is **conformal** to the hump (`η`-axis follows the wall) — handled automatically by `gridgen/benchmark/sweptwing_hump.m` |

The shipped defaults match the paper configuration; no edits needed to reproduce Fig. 10.

### Expected output

- Linear amplitude development (Fig. 10a) showing a slight stabilization downstream of the hump.
- Growth rate `α_i(x)` (Fig. 10c) with the characteristic dip-peak-dip signature around `x ∈ [820, 1000]`, matching the AHLNS reference.

### Runtime

- 2000×100 grid: ~0.6 h on a desktop (Intel Xeon, 4 cores, 32 GB)

---

## Notes for replication

- **Reference benchmark data.** NPSE (Cases 1, 2) and AHLNS (Case 3) reference curves are shipped in `literature/` and auto-loaded by the unified callers via `ref_file`. To override or skip the overlay, edit the `ref_file = ...` line near the top of the post-processing block in `callers/benchmark/<case>.m`.
- **Buffer treatment.** All three cases use parabolisation (`Opt.buffer = 'para'`) at the outflow. The buffer starts at 85% of the streamwise domain by default.
- **Convergence criterion.** `Opt.Conv = 1e-4` (default); for the nonlinear CFI case the absolute floor `1e-14` in `check_convergence.m` is what allows the harmonics to fully converge despite their tiny amplitudes.
- **Memory.** Cases 1 and 3 fit comfortably on a laptop. Case 2 (nonlinear CFI, 6 modes) needs roughly 15 GB of RAM at peak.
