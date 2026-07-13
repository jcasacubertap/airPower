# AmplitudeFinisher

Find the disturbance amplitude (`Stab.A0_fund`) that makes the **DeHNSSo
Harmonic Navier–Stokes** solution's spanwise-velocity perturbation match the
experimental **PIV**, for an OpenFOAM **flat-plate** base flow.

Pure MATLAB, driving DeHNSSo's `main()` directly. Two files matter:
**`AmplitudeFinisher.m`** (the caller you run) and
**`inputsAmplitudeFinisher.m`** (everything you set). Nothing else to touch.

## Run it

```matlab
cd .../instAbility/AmplitudeFinisher
AmplitudeFinisher
```

For a new case, edit the **MAIN INPUTS** block at the top of
`inputsAmplitudeFinisher.m`, then run `AmplitudeFinisher.m`. It writes
`AmplitudeFinisher_match.png` and `AmplitudeFinisher_match.mat` next to itself.

The whole airPower tree must be present (base-flow CSV, `DeHNSSo_input.mat`,
`inputs.jl`, PIV `Validation/`, and `airfoilFlowData/`). Sleep on the dev laptop
is handled by the global `local.keepawake` LaunchAgent.

## What it does

```
base-flow CSV ──gridgen──► StabGrid ──(truncate to the x/c window)──┐
PIV Gen/Case  ──────────►  w' target (raw plane, spanwise mean removed)
                                                                    │
                     findAmplitude:  linear seed → nonlinear sweep,
                     matched on the peak of the w'-RMS profile
                                                                    │
                                                                    ▼
                                        matched Stab.A0_fund + collapse plot
```

## MAIN inputs (top of `inputsAmplitudeFinisher.m`)
1. `gridgen.baseFlowCsv` — OpenFOAM mid-plane CSV (x,y,z,u,v,w,p,omz)
2. `Stab.lambda_z` — spanwise fundamental wavelength (7.5 mm)
3. `match.xcWindow` — chordwise stations to match, `[x/c_min x/c_max]` in %.
   Two close values (e.g. `[17 18]`) → a tight amplitude match. A wide range
   (e.g. `[17 30]`) → a **multi-station collapse test**: one A0 is fit over all
   those stations and every profile is plotted (see below).
4. `Stab.N` — spanwise Fourier modes (5 near the LE; **larger for the saturated
   tail** — the required N grows with amplitude)
5. `match.sweep` — amplitude bracket multipliers on the linear seed
6. `VAL.Gen/Case` — PIV case (`[]` → read from `airPower/inputs.jl`)
7. `machine` — **`'small'` | `'large'`** (see next section)

Everything below the MAIN block (reference scales, gridgen params, solver `Opt`,
PIV field names, frame `xOffset`) has sensible defaults.

## Small vs large machine (input #7)

`in.machine` sets the solver's memory/speed strategy in one place:

| `in.machine` | `lu_mode` | Newton takeover | use for |
|--------------|-----------|-----------------|---------|
| `'small'` | `'none'` (re-factorise each iter; low RAM, ~100 s/iter) | off | ≤16 GB laptop |
| `'large'` | `'full'` (cache all mode factorisations; ~10× faster, tens of GB) | on (deep-saturation convergence) | big-RAM workstation |

So on a big machine: set `in.machine = 'large';` and nothing else changes.
`'full'` at high N over the full domain needs tens of GB — that's what it's for.

## Multi-station collapse test

Set a **wide** `match.xcWindow`, e.g. `[17 30]`, and (on a large machine) raise
`Stab.N` and set `machine='large'`. `AmplitudeFinisher.m` then fits the single
A0 that best matches all those stations and plots each one. If HNS reproduces
the CFI growth, the profiles **collapse** onto the PIV at one amplitude; if not,
the per-station peak ratios drift (that drift = HNS-vs-experiment growth error).

Notes learned:
- Match on a **clean x/c window** — below ~16% the PIV has a near-wall artifact
  (its global-max sits at the wall, not the CFI peak); the deeply-saturated tail
  (x/c ≳ 27) needs **larger N** and is memory-heavy (→ `machine='large'`).
- **Streamwise frame**: OpenFOAM x is from the inlet, PIV x/c from the LE; they
  differ by `xInlet`. This is handled by `match.xOffset = xInlet` (airPower DFP
  convention `S = xInlet + x_OF`). Getting this wrong samples the base flow at
  the wrong station → a too-broad eigenfunction. Don't change it.
- The OpenFOAM base flow is the **true unperturbed laminar** flow (no MFD); PIV's
  spanwise mean carries the MFD downstream. The w' comparison removes the
  spanwise mean on both sides, so only the k≥1 perturbation is compared.
- Solver: `lu_mode='none'` fits ~16 GB but is slow; `'full'` is fast but needs
  big RAM (set via `machine`). Runs should be single-process (concurrent MATLABs
  caused memory-pressure kills on the laptop).

## Files
| file | role |
|------|------|
| `AmplitudeFinisher.m` | the caller — run this |
| `inputsAmplitudeFinisher.m` | MAIN + SECONDARY inputs |
| `src/baseflow/buildStabGrid.m` | OpenFOAM CSV → StabGrid (DeHNSSo gridgen, cached) |
| `src/baseflow/truncateStabGrid.m` | shorten domain to the matched window |
| `src/hns/runHNS.m` | one `main()` solve (linear/nonlinear) |
| `src/hns/extractWprimeHNS.m` | StabRes → w'-RMS profiles (base flow & MFD removed) |
| `src/validation/loadPIVwPrime.m` | PIV raw-plane w'-RMS target |
| `src/validation/stationArclength.m` | PIV x/c (% chord) → arc-length S |
| `src/validation/resolveGenCase.m` | read Gen/Case from `airPower/inputs.jl` |
| `src/matching/wPrimeObjective.m` | HNS-vs-PIV peak-RMS misfit (diagnostic) |
| `src/matching/findAmplitude.m` | linear seed → nonlinear peak-match sweep |

## Result so far (Gen0/Case0)
Frame-corrected peak match at x/c 17–18: **A0_fund ≈ 2.26×10⁻³** (peak + shape
match PIV; near-wall shoulder is PIV noise). A single-A0 collapse over the full
downstream range is the open item — needs a large machine (high N, `lu_mode`
`'full'`) to test whether HNS reproduces the CFI growth all the way.
