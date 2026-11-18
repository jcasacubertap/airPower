# airPower

Swept-wing crossflow-instability pipeline in three stages, driven by one
top-level dispatcher (`run.jl`) and one configuration, split over `inputs.jl`
(per case/run) and `settings.jl` (set once, rarely touched).

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

# instAbility — DeHNSSo stability solver        (case: cfi_flat_plate, cfi_hump, …)
julia run.jl instAbility DeHNSSo mesh cfi_flat_plate   # midPlane -> gridgen -> StabGrid -> grid figure
julia run.jl instAbility DeHNSSo run  cfi_flat_plate   # solver   -> StabRes (-> PostProcessing/io/input)

# PostProcessing — analysis                     (task: importData | reynoldsOrrProdTerms)
julia run.jl PostProcessing reynoldsOrrProdTerms
```

(`<case>` is the name of a DeHNSSo `gridgen/benchmark/<case>.m` + `callers/benchmark/<case>.m`
pair. Only `cfi_flat_plate` has its base flow staged from PreProcessing; other cases mesh the
base flow already in `baseflow/output/benchmark/`. `mesh` saves a grid-assessment figure to
`instAbility/io/plotting/DeHNSSo/<case>_grid.png`. `AmplitudeFinisher` is not in the dispatcher.)

## Data flow (coupling)

```
PreProcessing  postProcessing/midPlane.csv|.bin
      │  julia run.jl instAbility DeHNSSo mesh <case>
      ▼
   bf_<case>.csv ── gridgen ──► DeHNSSo_input_<case>.mat     (StabGrid = base flow)
      │  julia run.jl instAbility DeHNSSo run <case>
      ▼
   PostProcessing/io/input/<case>_<lin|nl>_N<N>_lz<λz>mm_A<A0>.mat   (StabRes + StabGrid)
      │  julia run.jl PostProcessing <task>
      ▼
   PostProcessing/io/output/<stab>_<task>.mat   (sRO + sBF + sPert + inp)
   PostProcessing/io/plotting/*.png
```

Filenames use `d` for `.` and `m` for `-` (e.g. `lz7.5mm` → `lz7d5mm`, `A5e-05` → `A5em05`).

## Configuration — `inputs.jl` + `settings.jl`

Two Julia files at the root, with the same block layout, form the single source
of truth:
- `inputs.jl` — what you set per case or per run (flow, geometry, resolution,
  bump, validation Gen/Case, PostProcessing load choices);
- `settings.jl` — set once and rarely touched (meshing details, export options,
  BL-metrics method, IBL validation grid, plot fine-tuning).

`PreProcessing/src/auxiliary/loadInputs.jl` merges them into one `inp`, so code
reads `inp.<Case>.<variable>` whichever file holds the entry. A key lives in
exactly one file; the same entry in both is an error. To carry a case's values
onto a newer airPower, use
`julia PreProcessing/src/selfRunning/seedInputs.jl <originCase> <destinationCase>`
(it handles both files, and an older single-file `inputs.jl` as origin).

PreProcessing reads them natively; for the MATLAB stages the dispatcher **exports** the relevant block to a
generated MATLAB config (`PostProcessing/inputs_gen.m`, gitignored) — MATLAB never
parses Julia. (DeHNSSo's own solver parameters currently live in its caller,
`instAbility/solvers/DeHNSSo/callers/benchmark/<case>.m`.)

Key `inp.PostProcessing` options:
- `task` — `importData` | `reynoldsOrrProdTerms` (the CLI `<task>` overrides it)
- `loadMode`, `fieldsFile` — the Stab input read from `io/input`
- `validation` — w-profile figures in x/c against PIV (`true`) or in S*/δ₀, DeHNSSo alone (`false`)
- in `settings.jl`: `ro.loadAnalysis` — `false`: compute + save + plot; `true`: load an `io/output` bundle and re-plot (no recompute); `plot.*` — the shared plot window

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
  inputs.jl         configuration: per case / run
  settings.jl       configuration: set once, rarely touched (merged with inputs.jl)
  PreProcessing/
    run.jl          Julia base-flow orchestrator (module dispatch)
    src/ modules/ …
  instAbility/
    solvers/DeHNSSo/          stability solver, v2 (gridgen + callers + src + docs)
    solvers/TrOY/             transient-growth / optimal-decay tooling
    applications/AmplitudeFinisher/   (standalone; not in the dispatcher)
    src/SolverCoupling/       airPower glue (grid-assessment figure)
    io/plotting/DeHNSSo/      <case>_grid.png                [gitignored]
  PostProcessing/
    run.m           MATLAB entry (loads the generated inputs_gen.m)
    src/            importData, reynoldsOrrProdTerms, plotters
    io/input/       Stab files (from `instAbility DeHNSSo run`)   [.mat gitignored]
    io/output/      analysis bundles {sRO, sBF, sPert, inp}       [.mat gitignored]
    io/plotting/    figures
```
