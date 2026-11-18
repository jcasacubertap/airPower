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

For other supported cases (`sweptwing_flat`, `sweptwing_hump`, `m3j` / OpenFOAM): see [README.md](../README.md).

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

### LU factorisation cache

The LHS matrix is factorised once per mode and reused across iterations. Three caching strategies:

```matlab
Opt.lu_mode = 'full';  % cache all LU factorisations (fastest, most memory)
Opt.lu_mode = 'auto';  % cache up to lu_max_cache, evict oldest (bounded memory)
Opt.lu_mode = 'none';  % solve with backslash each time (slowest, least memory)
% Opt.lu_max_cache = 3;  % only relevant for 'auto' mode
```

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
