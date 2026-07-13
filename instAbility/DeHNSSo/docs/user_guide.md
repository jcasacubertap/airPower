# DeHNSSo — User Guide

**Version 1.2** | Delft Harmonic Navier-Stokes Solver

---

## 1. What DeHNSSo does

DeHNSSo solves the nonlinear stability problem for laminar boundary layers. It computes the amplitude evolution of wave-like disturbances (Tollmien-Schlichting waves, cross-flow vortices) as they propagate downstream and interact nonlinearly.

The equations are the incompressible Navier-Stokes equations, written for perturbations on top of a steady base flow. The perturbations are expanded in a Fourier series in time and spanwise direction (harmonic balancing), and the resulting coupled system is solved iteratively until convergence.

---

## 2. Prerequisites

- MATLAB R2020b or later (no additional toolboxes required)
- A base-flow input file `DeHNSSo_input.mat` containing the `StabGrid` struct

---

## 3. Quick start

### Step 1 — Generate base flow and stability grid

For a Blasius flat-plate boundary layer:

```matlab
cd baseflow/IBL/
caller        % produces bf_blasius.mat in baseflow/output/
```

Then build the stability grid via the per-case script in `gridgen/benchmark/`:

```matlab
cd gridgen/benchmark/
blasius       % reads bf_blasius.mat from baseflow/output/benchmark/, writes DeHNSSo_input.mat
```

For other supported cases (`blasius_csv`, `sweptwing_flat`, `sweptwing_hump`, `m3j` / OpenFOAM): see [README.md](../README.md).

### Step 1b — Optional streamwise clustering (curvilinear cases)

By default the gridgen scripts distribute `nx` stations uniformly along the wall arc length. For cases with a localized geometric feature (hump, leading-edge curvature, roughness element) you can cluster stations around one or more centres via a Gaussian step-size modulation. Configured in the per-case gridgen script via three `params` fields:

| Parameter | Meaning |
|---|---|
| `params.mug` | Hump centre(s) in physical x (scalar or vector for multiple humps) |
| `params.sig` | Gaussian width as a fraction of the half-domain (default 0.2). Scalar broadcasts to all humps; vector must match `mug` length |
| `params.ag` | Step-size modulation at the peak: `<1` denser, `=1` no refinement, `>1` sparser (default 1). Scalar or vector |

Algorithm (consumed by `gridgen/src/resample_curvilinear.m`):

1. Map each `mug` into normalized ξ ∈ [-1, 1].
2. At every candidate ξ, compute `G_k(ξ) = exp(−½·((ξ − μ_n)/σ)²)` for each hump `k`.
3. Local modulation: `fg(ξ) = 1 − (1 − ag) · max_k G_k(ξ)` → `fg = ag` at each peak.
4. Build the new ξ grid by integrating `∫ dξ / fg` and renormalizing — denser where `fg` is small.

Example (heritage `sweptwing_hump.m`):

```matlab
x_m = 0.1837410879 / 2.1394e-4;   % physical hump centre / lref ≈ 859
params.mug = [x_m];               % single hump
params.sig = [0.2];               % Gaussian half-width fraction
params.ag  = [0.5];               % step size halved at the peak
```

Multi-hump example:

```matlab
params.mug = [x_m1, x_m2];        % two humps
params.sig = [0.2, 0.15];         % per-hump widths
params.ag  = [0.5, 0.7];          % per-hump strengths
```

Leave `params.mug = []` or `params.ag = 1` everywhere to fall back to uniform wall-arc-length spacing.

### Step 1c — Cartesian / OpenFOAM CSV input

The cartesian path (`gridgen/src/resample_cartesian.m`) reads either a heritage `.mat` (e.g. `bf_blasius.mat`) or a flat-lattice CSV from a CFD code (e.g. OpenFOAM `midPlane_*.csv`). For CSV input a few extra `params` fields are usually needed in the gridgen script:

| Parameter | Meaning |
|---|---|
| `params.rescale` | `true` if the CSV is dimensional (m, m/s). When set, X/Y are divided by `lref` and U/V/W by `Uref` before resampling |
| `params.Uref`, `params.lref` | Required when `rescale = true`. Together with kinematic viscosity ν they fix `params.Re = Uref·lref/ν` |
| `params.interp_method` | `interp2` method used on the rescaled grid. Default `'makima'` — Modified Akima, handles non-uniform wall-clustered CFD spacing without overshoot. Alternatives: `'spline'` (smoother interior but boundary ringing), `'linear'` (safest, lowest order), `'cubic'` (silently falls back to spline on non-uniform grids) |
| `params.xi_trim_inflow` | Fraction (0–1) of the x-span to drop at the inflow before resampling. Useful for CFD sources with noisy inflow BCs |
| `params.xi_trim_outflow` | Same at the outflow |

The path detects whether the CSV is structured (`X`, `Y` are 2-D arrays) or a flat lattice scatter (`X`, `Y` are paired column vectors of length `nx·ny`). The reshape from flat scatter to canonical `(ny × nx)` form is automatic — no manual preprocessing needed.

Example (heritage `blasius_csv.m`):

```matlab
params.rescale      = true;
params.Uref         = 20.9;                    % m/s
params.lref         = 2.688803650654011e-04;   % m
params.Re           = 3.719126161394363e+02;   % Uref·lref/ν_air
% params.interp_method   = 'makima';   % default; uncomment to override
% params.xi_trim_inflow  = 0.02;       % drop first 2% of x
% params.xi_trim_outflow = 0.0;
```

Leaving `params.n_eta_new` / `params.n_xi_new` empty keeps the source resolution. Setting them resamples to the requested size.

### Step 2 — Configure and run

Per-case stability scripts live in `callers/benchmark/` (`blasius.m`, `sweptwing_flat.m`, `sweptwing_hump.m`, `m3j.m`). All four share an identical sectioned layout — only case-specific values differ:

| Section | What to set |
|---|---|
| **Spectral content** | `Stab.M`, `Stab.N`, `Stab.omega_0`, `Stab.beta_0` |
| **Inflow conditions** | `Stab.IC`, `Stab.phaseRef`, `Opt.linear`, `Stab.A0_fund` |
| **Boundary conditions** | `Opt.bc_top`; optional `Stab.bcw`/`Stab.bct` for actuation/forcing |
| **Buffer** | `Opt.buffer` (`'on'`/`'para'`/`'off'`), `Opt.xb` |
| **Solver options** | `Opt.plot`, `Opt.Conv`, `Opt.TH` |
| **Performance** | `Opt.lu_mode`, `Opt.parfor`, `Opt.gpu` |
| **Nonlinear solver choice** | `Opt.solver` (`'picard'`/`'newton'`) |
| **Run** | `Opt.rerun`, call to `main()` |
| **Post-processing** | `plot_amplitudes`, `plot_stability` |

Run a caller from the MATLAB command window:

```matlab
cd callers/benchmark/
blasius           % or sweptwing_flat / sweptwing_hump / m3j
```

---

## 4. Spectral content

```matlab
Stab.M = 5;            % temporal modes: indices -M ... +M (total 2M+1)
Stab.N = 0;            % spanwise modes: indices -N ... +N (total 2N+1)
Stab.omega_0 = 0.034;  % fundamental angular frequency (non-dimensional)
Stab.beta_0  = 0;      % fundamental spanwise wavenumber (0 = 2D)
```

Total number of modes: `nf = (2M+1)(2N+1)`. Each mode `(m,n)` has frequency `m*omega_0` and wavenumber `n*beta_0`.

### Single oblique mode

To track a single oblique wave without spanwise harmonics, set `N = 0` with `beta_0 > 0`:

```matlab
Stab.M = 1;
Stab.N = 0;
Stab.beta_0 = 0.1;  % oblique wavenumber
```

The solver assigns `betavec = beta_0 * sign(Mvec)` for correct conjugate symmetry.

To get the index of mode `(m,n)` in the solution arrays:

```matlab
j = ModeToModeNumber(m, n, Stab.M, Stab.N);
```

---

## 5. Inflow conditions

```matlab
Stab.IC       = 'ILST';   % local stability eigenvalue problem at inflow
Stab.phaseRef = 'pwall';  % 'pwall' = wall pressure reference | 'umax' = max u
```

### Setting initial amplitudes

The simplest way is to set the fundamental amplitude as a scalar:

```matlab
Stab.A0_fund = 0.00125*sqrt(2);  % fundamental mode amplitude
```

`init.m` automatically builds the full amplitude vector with conjugate pairs. For per-mode control, set `Stab.A0` directly:

```matlab
Stab.A0 = zeros(1, (2*Stab.N+1)*(2*Stab.M+1));
Stab.A0(ModeToModeNumber( 1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);
Stab.A0(ModeToModeNumber(-1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);
```

Always set both the mode and its conjugate (index `-m`, `-n`) to ensure real physical fields. The amplitude is the max streamwise velocity perturbation, normalized by the freestream velocity.

---

## 6. Boundary conditions

By default all BCs are homogeneous (zero perturbation at wall and freestream). You can prescribe inhomogeneous BCs to model actuation or external disturbances.

### Arrays

```matlab
Stab.bcx          % (1 x nx_bc)       x-locations of BC profile
Stab.bcw(k,ix,j)  % (3 x nx_bc x nf)  wall BC:       k=1 u, k=2 v, k=3 w
Stab.bct(k,ix,j)  % (3 x nx_bc x nf)  freestream BC: k=1 u, k=2 v, k=3 w
```

### Wall actuation (blowing/suction)

Model a wall-normal velocity slot (e.g., to control TS waves):

```matlab
j_fund = ModeToModeNumber(1, 0, Stab.M, Stab.N);
j_conj = ModeToModeNumber(-1, 0, Stab.M, Stab.N);
x_dom  = StabGrid.xun;
x0_act = x_dom(1) + 0.3*(x_dom(end) - x_dom(1));  % actuator centre
sig    = 0.01*(x_dom(end) - x_dom(1));              % actuator half-width
A_act  = 1e-4;
Stab.bcx             = x_dom;
Stab.bcw             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
Stab.bcw(2,:,j_fund) = A_act * exp(-(x_dom - x0_act).^2 / sig^2);
Stab.bcw(2,:,j_conj) = conj(Stab.bcw(2,:,j_fund));
```

`k=2` is the wall-normal velocity. Set `k=1` for streamwise, `k=3` for spanwise.

### Freestream forcing

Prescribe a perturbation at `y=H` (e.g., an acoustic wave or free-stream turbulence):

```matlab
Stab.bct(1,:,j_fund) = A_fst * ones(1, numel(x_dom));  % u at y=H
Stab.bct(1,:,j_conj) = conj(Stab.bct(1,:,j_fund));
```

**Note:** for the `(0,0)` mean-flow distortion mode, the `v` freestream row carries the transpiration boundary condition and is not overridden by `bct`.

### Freestream BC type

`Opt.bc_top` is a cell array specifying the BC type per velocity component {u, v, w}:

```matlab
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};   % default: all zero at y=H

% BC types:
%   'H_DR'  = homogeneous Dirichlet  (q = 0)
%   'IH_DR' = inhomogeneous Dirichlet (q = bct value, requires Stab.bct)
%   'H_NM'  = homogeneous Neumann    (dq/dy = 0)
```

Pressure has no BC at the freestream (determined by the equations). The MFD `v` component is always free at the freestream (displacement velocity).

**Receptivity example** — acoustic wave forcing v at freestream:
```matlab
Opt.bc_top = {'H_DR', 'IH_DR', 'H_DR'};  % u=0, v=forced, w=0
Stab.bct(2,:,j_fund) = 1e-4 * ones(1, numel(x_dom));
Stab.bct(2,:,j_conj) = conj(Stab.bct(2,:,j_fund));
```

For backward compatibility, scalar strings `'dirichlet'` and `'neumann'` are auto-converted to `{'H_DR','H_DR','H_DR'}` and `{'H_NM','H_NM','H_NM'}` respectively.

The ILST inflow computation uses the same BC type.

---

## 7. Linear vs nonlinear mode

```matlab
Opt.linear = 'on';   % linear: NLT off, single iteration, AF = 1
Opt.linear = 'off';  % nonlinear (default): iterative solve with NLT
```

In linear mode the solver converges in one iteration (equivalent to linear PSE). Use this for parametric sweeps, N-factor computations, or any analysis where nonlinear interactions are not needed.

---

## 8. Solver options

### Buffer zone

```matlab
Opt.buffer = 'on';    % amplitude damping (prevents end reflections)
Opt.buffer = 'para';  % parabolisation (zeroes streamwise diffusion, LPSE-like)
Opt.buffer = 'off';   % no buffer
Opt.xb     = 85;      % buffer start as % of domain length
```

### Convergence

```matlab
Opt.Conv  = 1e-4;   % convergence criterion (relative amplitude change)
Opt.ConvF = 100;    % relaxation: use Opt.Conv*Opt.ConvF until AF=1
Opt.AFg   = 1.1;    % amplitude factor growth rate
Opt.AMAX  = 0.1;    % maximum amplitude factor
Opt.TH    = 1e-11;  % nonlinear term activation threshold
```

### Picard acceleration & Newton fallback

The nonlinear convergence loop is Picard iteration with AF (amplitude factor) continuation. Two optional accelerators speed it up; a JFNK Newton fallback handles cases where Picard stalls.

**Anderson acceleration** (Walker & Ni, 2011). A least-squares fixed-point accelerator that combines the last `m` Picard iterates to produce a faster update.

```matlab
Opt.anderson = 0;   % off (default — plain Picard)
Opt.anderson = m;   % history depth m (typical: 3-6)
```

How it behaves:

- Maintains a rolling buffer of the last `m` `(q_old, residual)` pairs.
- At each iteration, solves a small `m × m` least-squares problem (Tikhonov-regularised so it degrades gracefully when residuals are nearly parallel near saturation) and forms an accelerated step from the difference history.
- Safety fallback: if the accelerated step diverges (magnitude > 1.5× plain Picard), the iteration keeps the plain Picard step instead.
- The history is reset automatically when `AF` changes or the set of active modes changes (the fixed point itself has moved).

Often turns the linear convergence of plain Picard into a superlinear one near saturation, halving the iteration count on M=5 nonlinear cases.

**Newton takeover** (`Opt.newton_takeover`). When Picard stalls at AF < 1 (typically near saturation), DeHNSSo switches to a JFNK Newton solver to push to AF = 1.

```matlab
Opt.newton_takeover = 'on';   % default — auto-switch when Picard stalls
Opt.newton_takeover = 'off';  % stay in Picard regardless
Opt.newton_tol      = 1e-10;  % Newton convergence tolerance on ‖F‖
Opt.newton_precond  = 'block'; % 'block' | 'full'
```

The stall heuristic triggers when `dal` stops decreasing for ≥ 5 iters, or when `picard_af_count ≥ 15` at one AF, or on amplitude blow-up. Newton then converges the current AF using the existing LU cache as preconditioner and ramps further.

### LU mode

```matlab
Opt.lu_mode = 'lapack_band';  % LAPACK zgbsv/dgbsv MEX — fastest single solve (~3x);
                              %   needs ~6-12 GB; requires one-time MEX build (see below).
Opt.lu_mode = 'full';         % cache every LU factorisation (fastest reuse, unbounded memory)
Opt.lu_mode = 'auto';         % LRU-cached LUs, bounded memory (default for nonlinear runs)
Opt.lu_mode = 'none';         % backslash each time, no caching (slowest, baseline)
% Opt.lu_max_cache = 3;       % only relevant for 'auto'
```

**Picking quickly:**

| Run                                  | Recommended                                |
| ------------------------------------ | ------------------------------------------ |
| Linear single-mode (TS)              | `lapack_band`                              |
| Nonlinear multi-mode, plenty of RAM  | `auto` (caches reusable factors)           |
| Memory-constrained                   | `auto` with small `lu_max_cache`           |
| Validation / debugging               | `none` (reference, no caching)             |

### Building the `lapack_band` MEX (one-time setup per machine)

`lu_mode = 'lapack_band'` calls a MEX wrapper around LAPACK's banded LU
(`zgbsv`/`dgbsv`). The MEX binary is platform-specific and is **not** shipped
with the repo — each user builds it once.

```matlab
>> run('src/hns/mexbuild_dgbsv.m')
```

Requirements:

- MATLAB R2018a or later (uses the `-R2018a` interleaved-complex API).
- A C compiler configured for MEX (`mex -setup` once if you've never used MEX).
  - macOS: Xcode Command Line Tools.
  - Linux: `gcc`.
  - Windows: MSVC / MinGW.

If `lu_mode = 'lapack_band'` is selected but the MEX is missing, `solve_modes`
errors with a message pointing to `mexbuild_dgbsv.m`. Fall back to
`lu_mode = 'auto'` if you cannot build (e.g. no compiler available).

### Parallel and GPU acceleration

Both require the MATLAB Parallel Computing Toolbox. If omitted, auto-detected at startup.

```matlab
Opt.parfor = 'off';  % 'on'  parfor over LHS stations and modes (helps at large ny or many modes)
Opt.gpu    = 'off';  % 'on'  GPU solve via gpuArray \ (beneficial for ny >> 80)
```

**`parfor`:** parallelises the LHS station loop (Phase 1) and the mode solve loop (Phase 3). Memory overhead grows with the number of parallel workers — each worker holds a copy of the broadcast LU objects. Effective crossover is roughly ny ≥ 80 where per-station compute exceeds broadcast cost.

**`gpu`:** transfers each assembled mode matrix `ML` to GPU memory and solves with `gpuArray \`. GPU does not support the `decomposition` cache; the matrix is refactorised on every iteration. Beneficial when ny is large (O(ny²) compute dominates memory bandwidth).

### Output control

```matlab
Opt.plot = 'on';   % 'on'  show amplitude and convergence plots during iteration (default)
Opt.plot = 'off';  % 'off' suppress all intermediate plots (recommended for batch runs)
```

---

## 9. Running and reruns

```matlab
Opt.rerun = 'off';  % full run: builds LHS matrices and returns them in LHS struct
Opt.rerun = 'on';   % reuse LHS matrices from a previous run (faster)
```

For a single run:

```matlab
[StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
```

For a parameter sweep (same base flow, different omega):

```matlab
omega_vec = Stab.omega_0 * linspace(0.8, 1.2, 5);
for k = 1:numel(omega_vec)
    Stab.omega_0 = omega_vec(k);
    if k == 1
        [StabRes_sweep{k}, ~, LHS] = main(Stab, StabGrid, Opt);
    else
        Opt.rerun = 'on';
        [StabRes_sweep{k}, ~, LHS] = main(Stab, StabGrid, Opt, LHS);
    end
end
```

When `Opt.rerun = 'on'`:
- LHS matrices (`MLA1`, `MLA2`, `MFD_diff`) are reused — no reassembly.
- LU factorisations are reused when `omega` and `beta` are unchanged (actuation/amplitude sweeps).
- LU factorisations are recomputed when `omega` or `beta` changes (frequency/wavenumber sweeps).

Indicative cost on the nonlinear Blasius TS case (M=5, Nx=800, Ny=40, `lu_mode='full'`):

| call | what runs | wall-clock |
|---|---|---|
| first (`rerun='off'`) | assembly + 6 LU factorisations + Picard | ~40 s |
| reuse (`rerun='on'`, amplitude sweep) | restored blocks + 6 cached factorisations, Picard only | ~6–8 s (~6× faster) |

The reused solve reproduces the from-scratch result to machine precision.

---

## 10. Output

The main output is `StabRes`:

| Field | Size | Description |
|---|---|---|
| `StabRes.A(j,i)` | `nf × nx` | Mode amplitude at each station |
| `StabRes.alpha(j,i)` | `nf × nx` | Streamwise wavenumber |
| `StabRes.u(j,y,i)` | `nf × ny × nx` | Streamwise velocity shape function |
| `StabRes.v(j,y,i)` | `nf × ny × nx` | Wall-normal velocity shape function |
| `StabRes.w(j,y,i)` | `nf × ny × nx` | Spanwise velocity shape function |
| `StabRes.p(j,y,i)` | `nf × ny × nx` | Pressure shape function |

Index convention: `y=1` is the freestream (`y=H`), `y=ny` is the wall (`y=0`).

### Figures produced during a run

| Figure | Content |
|---|---|
| 1 | Inflow: base flow profiles, eigenvalue spectrum, selected shape function |
| 2 | Amplitude evolution during convergence iterations |
| 3 | Convergence history (relative error per mode, log scale) |
| 4 | Final amplitude development (optionally overlaid with PSE/DNS reference) |
| 5 | N-factor, growth rates, wavenumber, phase speed |
| 6 | Eigenfunctions at selected streamwise stations |

### Post-processing

The unified callers in `callers/benchmark/` already wire post-processing for you:

```matlab
% Reference for overlay (set ref_file = [] to skip)
ref_file  = fullfile(rootdir, 'literature', 'StabRes_<case>_<METHOD>.mat');
ref_label = '<METHOD>';

if strcmpi(Opt.plot, 'on')
    if ~isempty(ref_file) && isfile(ref_file)
        plot_amplitudes(StabGrid, StabRes, struct('file', ref_file, 'label', ref_label));
    else
        plot_amplitudes(StabGrid, StabRes);
    end
    plot_stability(StabGrid, StabRes, Opt);
end
```

Reference data shipped in `literature/`:

| File | Used by | Source |
|---|---|---|
| `StabRes_Blasius_NPSE.mat` | `blasius.m` | NPSE reference (Westerbeek 2024 Fig. 8) |
| `StabRes_SweptWing_flat_NPSE.mat` | `sweptwing_flat.m` | NPSE reference (Westerbeek 2024 Fig. 9) |
| `StabRes_SweptWing_Hump_AHLNS.mat` | `sweptwing_hump.m` | AHLNS reference (Westerbeek 2024 Fig. 10) |

To overlay a different reference, point `ref_file` at any `.mat` file with a `StabRes` struct in the format that `plot_amplitudes` expects (see `src/plot/plot_amplitudes.m`).

---

## 11. Common issues

| Problem | Likely cause | Fix |
|---|---|---|
| ILST does not find a physical mode | Domain too short/tall, wrong `omega`/`beta` | Check eigenvalue spectrum in Figure 1 |
| Solver does not converge | Amplitude too high, too few modes | Reduce `Stab.A0`, increase `Stab.M` |
| Memory error | Large `nx × ny`, `lu_mode = 'full'` | Use `lu_mode = 'auto'` with `lu_max_cache`; disable `parfor` to avoid worker copies |
| Results differ from PSE | Buffer too short, nonlinear effects | Increase `Opt.xb`, check `Stab.M` |
| Freestream BC oscillations | Domain truncated, Dirichlet BC | Use `'H_NM'` on affected components: e.g. `Opt.bc_top = {'H_NM','H_NM','H_NM'}` |
