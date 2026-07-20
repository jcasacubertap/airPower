# airPower

Swept-wing crossflow-instability pipeline in three stages, driven by one
top-level dispatcher (`run.jl`) and one config file (`inputs.jl`).

| stage | language | role |
|-------|----------|------|
| `PreProcessing`  | Julia  | OpenFOAM base-flow generation (DirectFlatPlate / TunnelToCurvedPlate) |
| `instAbility`    | MATLAB | DeHNSSo harmonic Navier–Stokes stability solver (base flow → perturbation) |
| `PostProcessing` | MATLAB | analysis of the base flow + perturbation (Reynolds–Orr production, …) |

## Running — `julia run.jl <stage> …`

The top-level `run.jl` dispatches to each stage, shelling out to `matlab` for the
MATLAB stages. If `matlab` is not on your `PATH` (e.g. macOS), set the executable
via `AIRPOWER_MATLAB`, e.g. `AIRPOWER_MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab`.

```bash
# PreProcessing — OpenFOAM base flow            (action: all | clean | prep | viz | …)
julia run.jl PreProcessing DirectFlatPlate all
julia run.jl PreProcessing TunnelToCurvedPlate all

# instAbility — DeHNSSo stability solver        (case: sweptwing_flat)
julia run.jl instAbility DeHNSSo mesh sweptwing_flat   # midPlane -> gridgen -> StabGrid
julia run.jl instAbility DeHNSSo run  sweptwing_flat   # solver   -> StabRes (-> PostProcessing/io/input)

# PostProcessing — analysis                     (task: importData | reynoldsOrrProdTerms)
julia run.jl PostProcessing reynoldsOrrProdTerms
```

(`instAbility DeHNSSo` is wired for `sweptwing_flat`; `m3j` and `AmplitudeFinisher`
are not in the dispatcher yet.)

## Data flow (coupling)

```
PreProcessing  postProcessing/midPlane.csv|.bin
      │  julia run.jl instAbility DeHNSSo mesh <case>
      ▼
   bf_<case>.csv ── gridgen ──► DeHNSSo_input.mat            (StabGrid = base flow)
      │  julia run.jl instAbility DeHNSSo run <case>
      ▼
   PostProcessing/io/input/<case>_<lin|nl>_N<N>_lz<λz>mm_A<A0>.mat   (StabRes + StabGrid)
      │  julia run.jl PostProcessing <task>
      ▼
   PostProcessing/io/output/<stab>_<task>.mat   (sRO + sBF + sPert + inp)
   PostProcessing/io/plotting/*.png
```

Filenames use `d` for `.` and `m` for `-` (e.g. `lz7.5mm` → `lz7d5mm`, `A5e-05` → `A5em05`).

## Configuration — `inputs.jl`

One Julia file at the root is the single source of truth. PreProcessing reads it
natively; for the MATLAB stages the dispatcher **exports** the relevant block to a
generated MATLAB config (`PostProcessing/inputs_gen.m`, gitignored) — MATLAB never
parses Julia. (DeHNSSo's own solver parameters currently live in its caller,
`instAbility/DeHNSSo/callers/benchmark/<case>.m`.)

Key `inp.PostProcessing` options:
- `task` — `importData` | `reynoldsOrrProdTerms` (the CLI `<task>` overrides it)
- `loadMode`, `fieldsFile` — the Stab input read from `io/input`
- `loadAnalysis` — `false`: compute + save + plot; `true`: load an `io/output` bundle and re-plot (no recompute)
- `ro.*` — Reynolds–Orr options (modes, integration window, plot buffer fraction)

### PostProcessing: two entry points, both honour `inputs.jl`

- **Terminal / dispatcher** — `julia run.jl PostProcessing <task>`: uses the CLI
  `<task>`, saves the bundle + PNGs. The dispatcher writes `inputs_gen.m` and sets
  `AIRPOWER_PP_DISPATCH`, so `run.m` loads it as-is.
- **Directly in MATLAB** — open `PostProcessing/run.m` and Run it (or type `run`):
  leaves `sBF`/`sPert`/`sRO`/`inp` in the workspace for interactive work. With no
  dispatch flag set, `run.m` first re-syncs `inputs_gen.m` from `inputs.jl` (shells
  `julia … PostProcessing config`), so it always mirrors the current `inputs.jl`
  (task = its default). Julia is found via `AIRPOWER_JULIA` → `PATH` →
  `~/.juliaup/bin/julia`.

## Layout

```
airPower/
  run.jl            top-level dispatcher
  inputs.jl         single config source of truth
  PreProcessing/
    run.jl          Julia base-flow orchestrator (module dispatch)
    src/ modules/ …
  instAbility/
    DeHNSSo/        stability solver (gridgen + callers)
    AmplitudeFinisher/   (standalone; not yet in the dispatcher)
  PostProcessing/
    run.m           MATLAB entry (loads the generated inputs_gen.m)
    src/            importData, reynoldsOrrProdTerms, plotters
    io/input/       Stab files (from `instAbility DeHNSSo run`)   [.mat gitignored]
    io/output/      analysis bundles {sRO, sBF, sPert, inp}       [.mat gitignored]
    io/plotting/    figures
```
