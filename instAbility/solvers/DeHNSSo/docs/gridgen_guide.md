# gridgen — User Guide

gridgen turns a base flow into the `StabGrid` struct the solver reads: it resamples onto the
numerical grid, clusters the mesh where you need resolution, computes the base-flow
derivatives, and stores any wall forcing. It is a separate module (`gridgen/`) and runs
before the solver, once per case.

Per-case scripts live in `gridgen/benchmark/`. Run one, then run the matching caller in
`callers/benchmark/`.

## Contents

- [1. Generating a grid](#1-generating-a-grid)
- [2. Streamwise clustering](#2-streamwise-clustering)
- [3. Cartesian / OpenFOAM CSV input](#3-cartesian--openfoam-csv-input)
- [4. Base-flow derivatives](#4-base-flow-derivatives)
- [5. Wall forcing](#5-wall-forcing)
- [6. Steps and blocks on the wall](#6-steps-and-blocks-on-the-wall)

---

## 1. Generating a grid

For a Blasius flat-plate boundary layer:

```matlab
cd baseflow/caller/
caller_IBL    % writes bf_ts_blasius.mat to baseflow/output/ (scratch, gitignored)
```

Then build the stability grid via the per-case script in `gridgen/benchmark/`:

```matlab
cd gridgen/benchmark/
ts_blasius       % reads the COMMITTED bf_ts_blasius.mat in baseflow/output/benchmark/, writes DeHNSSo_input_ts_blasius.mat
```

The benchmarks read the committed base flows in `baseflow/output/benchmark/`, not the caller's
scratch output: to run a benchmark on a regenerated base flow, copy it there first. The other
cases (`ts_compressible`, `cfi_flat_plate`, `cfi_hump`, `cfi_step`, `cfi_curved_airfoil`, `resonator`,
`wall_jet`, `newton_continuation`) follow the same pattern: see [README.md](../README.md).

Every per-case script fills one `params` struct. The basic fields, all optional:

| Parameter | Meaning | Default |
|---|---|---|
| `params.n_xi_new`, `params.n_eta_new` | number of streamwise and wall-normal points | source resolution |
| `params.H` | domain height | `BF.H`, else 0.99 of the base-flow height |
| `params.y_i` | Malik mapping parameter: half of the wall-normal points lie below it | `BF.y_i`, else `H/10` |
| `params.xi_range` | `[lo, hi]` fractions of the streamwise span to keep | whole span |
| `params.FD_xi_order_1`, `params.FD_xi_order_2` | order of the streamwise first and second derivatives of the metrics | 4, 2 |
| `params.FD_eta_method`, `params.FD_eta_order` | wall-normal operator: `cheb`, or `fd` with this order | `'cheb'`, 4 |
| `params.snap` | move a node onto the box edges of a single step (see [6.](#6-steps-and-blocks-on-the-wall)) | `true` |
| `params.elliptic`, `params.elliptic_opts` | elliptic smoothing of an unstructured body-fitted mesh; `elliptic_opts` overrides `max_iter`, `tol`, `omega`, `n_uniform` | `false` |
| `params.plot`, `params.plot_opts` | show the raw and resampled base flow; `plot_opts.xi_slices` / `eta_slices` pick the profiles | `false` |

## 2. Streamwise clustering

By default the gridgen scripts distribute `nx` stations uniformly in ξ (wall arc length for curvilinear cases, `x` for the flat plate). For cases with a localized feature (hump, leading-edge curvature, roughness element) — or simply to concentrate stations where the disturbance grows fastest — you can cluster stations around one or more centers via a Gaussian step-size modulation. Configured in the per-case gridgen script via three `params` fields:

| Parameter | Meaning |
|---|---|
| `params.mug` | Clustering center(s) in physical x (scalar or vector for multiple centers). Taken from `params.forcing.xc` when a forcing is defined, so set it explicitly only for a feature that is not a wall source, such as a hump |
| `params.spread` | How far the refinement reaches: half-width in physical x. Defaults to the `params.forcing` half-width. Larger means gentler stretching |
| `params.refine` | How many times finer Δξ is at the peak than on an equidistant mesh |

All three are in physical units, so they read the same way as the feature they follow. Internally `spread` becomes the Gaussian width `sig = 3·spread / half-domain` and `refine` becomes `ag = 1/refine`; `params.sig` and `params.ag` are **not** accepted as inputs and raise an error naming the replacement.

These fields work for both the curvilinear (`resample_curvilinear.m`) and flat-plate (`resample_cartesian.m`) paths, which share one implementation in `gridgen/src/resample/cluster_xi.m`.

Algorithm (`gridgen/src/resample/cluster_xi.m`):

1. Map each `mug` into normalized ξ ∈ [-1, 1].
2. At every candidate ξ, compute `G_k(ξ) = exp(−½·((ξ − μ_n)/σ)²)` for each hump `k`.
3. Local modulation: `fg(ξ) = 1 − (1 − ag) · max_k G_k(ξ)` → `fg = ag` at each peak.
4. Build the new ξ grid by integrating `∫ dξ / fg` and renormalizing — denser where `fg` is small.

Example (`gridgen/benchmark/cfi_hump.m`):

```matlab
x_m = 0.1837410879 / 2.1394e-4;   % physical hump center / lref ≈ 859
params.mug    = [x_m];            % single hump (no forcing here, so set the centre)
params.spread = 55.137096;        % half-width the refinement reaches
params.refine = 2;                % dxi twice as fine at the peak
```

Multi-hump example:

```matlab
params.mug    = [x_m1, x_m2];     % two humps
params.spread = [55.1, 41.4];     % per-hump half-widths
params.refine = [2, 1.43];        % per-hump strengths
```

Gridgen reports what the mesh actually does, not just what was asked for: the finest Δξ as a multiple of the equidistant spacing, and the stretching, the worst ratio between neighbouring Δξ (1.00 means none). The measured refinement is always milder than the nominal `refine`, because compressing the mesh at the peak stretches it elsewhere. Aim for stretching near 1.0; widen `spread` or lower `refine` if it grows.

Leave `params.mug = []` and define no `params.forcing`, or set `params.refine = 1`, to fall back to uniform spacing (wall arc length for curvilinear, `x` for the flat plate).

## 3. Cartesian / OpenFOAM CSV input

The cartesian path (`gridgen/src/resample/resample_cartesian.m`) reads either a heritage `.mat` (e.g. `bf_ts_blasius.mat`) or a flat-lattice CSV from a CFD code (e.g. OpenFOAM `midPlane_*.csv`). For CSV input a few extra `params` fields are usually needed in the gridgen script:

| Parameter | Meaning |
|---|---|
| `params.rescale` | `true` if the CSV is dimensional (m, m/s). When set, X/Y are divided by `lref` and U/V/W by `Uref` before resampling |
| `params.Uref`, `params.lref` | Required when `rescale = true`. Together with kinematic viscosity ν they fix `params.Re = Uref·lref/ν` |
| `params.interp_method` | `interp2` method used on the rescaled grid. Default `'makima'` — Modified Akima, handles non-uniform wall-clustered CFD spacing without overshoot. Alternatives: `'spline'` (smoother interior but boundary ringing), `'linear'` (safest, lowest order), `'cubic'` (silently falls back to spline on non-uniform grids) |
| `params.xi_trim_inflow` | Fraction (0–1) of the x-span to drop at the inflow before resampling. Useful for CFD sources with noisy inflow BCs |
| `params.xi_trim_outflow` | Same at the outflow |

The path detects whether the CSV is structured (`X`, `Y` are 2-D arrays) or a flat lattice scatter (`X`, `Y` are paired column vectors of length `nx·ny`). The reshape from flat scatter to canonical `(ny × nx)` form is automatic — no manual preprocessing needed.

Example, a dimensional CSV at 20.9 m/s non-dimensionalised on the Blasius length at its inflow:

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

## 4. Base-flow derivatives

`gridgen` computes the physical base-flow derivatives that the stability operator needs and stores them on the grid, so the solver never differentiates the base flow itself:

- First derivatives `∂/∂x`, `∂/∂y` of `U, V, W` (and of `rho, T` for a compressible base flow).
- Compressible only: second streamwise `∂²/∂x²` (`dxxU/V/W/T`), second wall-normal `∂²/∂y²` (`ddyU/V/W/T`), and mixed `∂²/∂x∂y` (`dxyU`, `dxyV`).

By default these come from finite differences of the (interpolated) base flow. If you have an **analytic** base flow you can instead supply the derivatives directly on the BF struct — gridgen detects them automatically and skips the corresponding FD differentiation. No flag is needed:

```matlab
BF.dxU = ...;  BF.dyU = ...;   % supply both halves of a pair
```

Rules:

- Detection is **per field** and **per pair**: supply *both* `dx<f>` and `dy<f>` to replace that field's first-derivative differencing. Supply only one and gridgen warns (`<f>: dx<f> provided without its pair ...`) and recomputes both by FD.
- Recognized fields: `dxU,dyU, dxV,dyV, dxW,dyW` (incompressible); plus `dxrho,dyrho, dxT,dyT` and the second derivatives `dxxU/V/W/T, ddyU/V/W/T, dxyU,dxyV` (compressible). The single second-derivative fields have no pair, so they are taken as-is when present.
- gridgen reports what it used: `compute_bf_derivatives: using N provided base-flow derivative field(s)`.

Supplying analytic derivatives avoids the finite-difference noise that interpolated CFD base flows can introduce near the wall (which otherwise inflates the wall shear by `O(N²)`). For a smooth analytic base flow the FD path is already accurate to ~0.1%, so this is mainly useful for noisy or analytically-known base flows.


---

## 5. Wall forcing

A wall source is declared here, in the gridgen script: gridgen refines the mesh to resolve
it and stores the profile on the grid. The caller then sets only its strength. Location,
width and shape are therefore defined once. For the caller side, see `docs/DeHNSSo_guide.md`,
"Wall actuation".

`params.forcing` describes what sits at the wall. Two choices decide the form: where the **shape** comes from (the built-in cosine, or an array you supply), and **how many** sources there are (one, or an array of them). That gives seven forms:

| # | shape | how many | spec |
|---|---|---|---|
| 1 | cosine | one | `struct('xc', xc, 'w', w)` |
| 2 | cosine | N, same width | `struct('xc', [xc1 xc2 xc3], 'w', w)` |
| 3 | cosine | N, own widths | `struct('xc', [xc1 xc2], 'w', [w1 w2])` |
| 4 | tabulated | one | `struct('xc', xc, 'w', w, 'v', v_data)` |
| 5 | tabulated | one, uneven samples | `struct('xc', xc, 'w', w, 'x', x_data, 'v', v_data)` |
| 6 | tabulated | N | `struct('xc', [xc1 xc2 xc3], 'w', w, 'v', v_data)` |
| 7 | tabulated | one, **placed by the data** | `struct('x', x_data, 'v', v_data)` |

`'shape','cosine'` may be given explicitly in 1-3, but it is the default and `'cosine'` is the only value. In 4-6 `v_data` may be complex, carrying phase.

**Forms 1-6 place the source themselves.** `xc` is the centre and `w` the HALF-width, so the footprint spans `2w`; a supplied array is stretched onto `[xc-w, xc+w]` and normalised to unit peak. Its own units, origin and sample order therefore never enter, and in form 5 `x_data` is read for the *relative spacing* of the samples alone. Pass a profile tabulated in millimetres about its own centre exactly as it comes out of the file.

**Form 7 is the exception**: with no `xc`/`w`, the centre and half-width are taken from where `|v|` first and last exceeds `1e-6` of its peak, so there `x_data` **is** the position and must already be in the grid's coordinate. Use it when the distribution is the specification and has no natural centre and width to quote. Two consequences: a profile that is exactly zero at its end points reports a slightly narrow `w`, and a units mistake moves the source rather than being ignored.

Whichever form, a source lying entirely off the domain raises `main_gridgen:forcingOffDomain`, and one hanging over an edge warns with `forcingClipped`.

```matlab
% one actuator
params.forcing = struct('xc', 1379.1, 'w', 18);

% a measured neck profile, as it comes out of the file
D = load('hAbs_xProfile_TM295hz.mat');           % h is a shape; its scale is discarded
params.forcing = struct('xc', Xc, 'w', s_nd/2, 'v', D.h);

% an array of three identical necks
params.forcing = struct('xc', Xc + (0:2)*(s_nd + d_nd), 'w', s_nd/2, 'v', D.h);

params.spread = 30;     % how far the refinement reaches
params.refine = 5;      % dxi at the centre vs equidistant
```

`spread` and `refine` are the same parameters as in [2. Streamwise clustering](#2-streamwise-clustering). `spread` defaults to the forcing half-width, which is right for an actuator; set it explicitly only when the refinement must reach further than the source itself, as a resonator neck does.

**Sources that differ from each other** go in a cell array, one struct each. This composes with the table: a vector `xc` inside a cell entry still expands, so the example below is four sources.

```matlab
params.forcing = { struct('xc', [xc1 xc2 xc3], 'w', 18, 'v', v_data), ...   % three alike
                   struct('xc', xc4, 'w', 24) };                            % one different
```

gridgen stores `StabGrid.forcing` (the unit-peak profile(s) sampled on the refined grid) plus `StabGrid.cluster.xc/.w`. The distribution's `v` values need only encode the **shape**; the amplitude is set in the caller. Note the fixed station budget is shared across all centers, so raise `params.n_xi_new` when refining several spots.

---

## 6. Steps and blocks on the wall

A step or a block on the wall is not meshed. The wall stays flat and gridgen stores the solid as one box `[x_lo, x_hi, h]` per step. The solver masks the boxes out. For the solver side, see `docs/DeHNSSo_guide.md`, "Steps and blocks on the wall".

gridgen finds the boxes in the base flow: a solid computed on a Cartesian grid shows up as nodes with zero velocity above the wall. Nothing about the geometry is typed into the script.

```
Solid detected: STEPS, 2 box(es), 2 face(s)
   x =  79.8750 .. 149.8750,  h = 1.2000
   x = 149.8750 .. 219.8750,  h = 2.4000
```

| Parameter | Meaning |
|---|---|
| `params.y_i` | Set it to the step height, so about half the wall-normal points lie below the step |
| `params.mug`, `params.refine`, `params.spread` | Streamwise refinement at the face, as in [2. Streamwise clustering](#2-streamwise-clustering). Use `params.mug = StabGrid.step.x` for the first face |
| `params.StepType`, `params.StepX`, `params.StepH` | Optional. State a single forward-facing step by hand, all three together. gridgen still detects the step and warns if the two disagree |

What gridgen writes:

| Field | Meaning |
|---|---|
| `StabGrid.step.type` | `FFS`, `BFS`, `BLOCK` or `STEPS` |
| `StabGrid.step.seg` | the boxes, one row `[x_lo, x_hi, h]` each |
| `StabGrid.step.seg_tol` | how well each box is known, `[x_tol, h_tol]` |
| `StabGrid.step.snapped` | true when the grid was put exactly on the box |

Limits to keep in mind:

- Only Cartesian input is read. On a curvilinear or unstructured base flow a solid would be treated as fluid, so gridgen warns if it sees one.
- A face lies somewhere between the last fluid column and the first solid one. gridgen places it halfway and reports half that gap as `x_tol`. The height is read off the base-flow nodes, so it is exact only when the surface falls on a node.
- A solid that does not rest on the wall, or one that spans the whole domain, is refused.
