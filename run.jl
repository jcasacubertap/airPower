#!/usr/bin/env julia
#
# airPower — top-level dispatcher
#
# Routes a stage command to the right runner. The stages are mixed-language:
# PreProcessing is Julia; instAbility and PostProcessing are MATLAB (launched via
# `matlab -batch` under the hood). This dispatcher shells out to each stage so
# they stay isolated in their own process.
#
# Usage:
#   julia run.jl PreProcessing <module> <action>          # e.g. DirectFlatPlate all
#   julia run.jl instAbility   DeHNSSo <mesh|run> <case>  # (wired in a later phase)
#   julia run.jl PostProcessing <task>                    # importData | reynoldsOrrProdTerms
#
# MATLAB executable: resolved automatically — $AIRPOWER_MATLAB, else `matlab` on
# PATH, else the newest standard install (/Applications/MATLAB_R*.app on macOS,
# /usr/local/MATLAB/R* on Linux). Set AIRPOWER_MATLAB to force a particular one.
#

const ROOT = @__DIR__
include(joinpath(ROOT, "inputs.jl"))   # `inp` — the single source of truth

function usage()
    println("""
airPower dispatcher

  julia run.jl PreProcessing <module> <action>
      module  : DirectFlatPlate | TunnelToCurvedPlate
      action  : all | clean | prep | viz | ...
  julia run.jl instAbility DeHNSSo <mesh|run> <case>
      mesh    : PreProcessing midPlane -> gridgen -> StabGrid
      run     : DeHNSSo caller -> StabRes (-> PostProcessing/io/input/<encoded>.mat)
      case    : sweptwing_flat
  julia run.jl PostProcessing <task>
      task    : importData | reynoldsOrrProdTerms
""")
end

# --- MATLAB launcher --------------------------------------------------------
# Resolution ladder, so a fresh clone runs with no environment setup at all:
#   1. $AIRPOWER_MATLAB       explicit override, wins outright
#   2. `matlab` on PATH       the usual properly-linked install
#   3. standard install dirs  newest release first
# Step 3 is not redundant. On macOS `matlab` is very often a SHELL ALIAS rather
# than a real executable, and an alias is invisible to a spawned process: the
# dispatcher would die with ENOENT even though typing `matlab` works fine.
const MATLAB_EXE = Ref("")

# Where a standard MATLAB install lives, per platform. The layout is the same
# everywhere — <root>/<release dir>/bin/matlab[.exe] — only the roots differ.
function matlab_roots()
    if Sys.isapple()
        ["/Applications"]                                   # MATLAB_R2025b.app
    elseif Sys.iswindows()
        pf   = get(ENV, "PROGRAMFILES",      "C:\\Program Files")
        pf86 = get(ENV, "PROGRAMFILES(X86)", "C:\\Program Files (x86)")
        [joinpath(pf, "MATLAB"), joinpath(pf86, "MATLAB")]  # ...\MATLAB\R2025b
    else
        ["/usr/local/MATLAB", "/opt/MATLAB", "/opt/matlab"] # .../R2025b
    end
end

matlab_exename() = Sys.iswindows() ? "matlab.exe" : "matlab"

function find_matlab()
    isempty(MATLAB_EXE[]) || return MATLAB_EXE[]

    override = get(ENV, "AIRPOWER_MATLAB", "")
    if !isempty(override)
        isfile(override) || error("AIRPOWER_MATLAB is set to '$override', which is not a file.")
        return MATLAB_EXE[] = override
    end

    onpath = Sys.which("matlab")
    onpath === nothing || return MATLAB_EXE[] = onpath

    roots = matlab_roots()
    found = Tuple{String,String}[]                    # (release tag, executable)
    for r in roots
        isdir(r) || continue
        for d in (try readdir(r) catch; String[] end)
            exe = joinpath(r, d, "bin", matlab_exename())
            isfile(exe) || continue
            tag = match(r"R(\d{4}[ab])", d)
            push!(found, (tag === nothing ? "" : String(tag.captures[1]), exe))
        end
    end
    isempty(found) && error("""
        MATLAB not found. Tried, in order:
          1. \$AIRPOWER_MATLAB  — unset
          2. `matlab` on PATH    — not found. NOTE a shell alias does not count:
                                   it is invisible to a spawned process. A MATLAB
                                   launched from Finder/Dock/Start menu also sees
                                   only a minimal environment.
          3. $(join(roots, ", "))  — nothing matching */bin/matlab
        Point AIRPOWER_MATLAB at the binary, e.g.
          export AIRPOWER_MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab
        """)

    sort!(found; by = first, rev = true)               # newest release first
    exe = found[1][2]
    @info "airPower: using MATLAB at $exe (set AIRPOWER_MATLAB to override)"
    return MATLAB_EXE[] = exe
end

# Extra environment for the child is passed as keywords, e.g.
#   matlab_batch(stmt; AIRPOWER_PP_DISPATCH = "1")
function matlab_batch(statement; env...)
    cmd = `$(find_matlab()) -batch $statement`
    isempty(env) && return run(cmd)
    run(addenv(cmd, [String(k) => v for (k, v) in env]...))
end

# --- serialize a Julia value to a MATLAB literal ---
mat_lit(x::AbstractString) = "'" * x * "'"
mat_lit(x::Bool)           = x ? "true" : "false"
mat_lit(x::Real)           = string(x)
mat_lit(x::AbstractVector) = isempty(x) ? "[]" : "[" * join(string.(x), " ") * "]"

# --- PreProcessing: shell to the Julia orchestrator in its own process ---
function run_preprocessing(args)
    script = joinpath(ROOT, "PreProcessing", "run.jl")
    @info "airPower ▶ PreProcessing" args = args
    run(`$(Base.julia_cmd()) $script $args`)
end

# --- instAbility / DeHNSSo: mesh (gridgen) and run (solver) ---
function run_instability(args)
    length(args) >= 3 ||
        error("usage: julia run.jl instAbility DeHNSSo <mesh|run> <case>")
    modul, sub, case_ = args[1], args[2], args[3]
    modul == "DeHNSSo" || error("instAbility: only 'DeHNSSo' is wired (got '$modul').")
    dehnsso = joinpath(ROOT, "instAbility", "DeHNSSo")
    if sub == "mesh"
        dehnsso_mesh(dehnsso, case_)
    elseif sub == "run"
        dehnsso_run(dehnsso, case_)
    else
        error("instAbility DeHNSSo: subcommand must be 'mesh' or 'run' (got '$sub').")
    end
end

# mesh: copy PreProcessing midPlane -> base-flow csv, then run gridgen -> StabGrid
function dehnsso_mesh(dehnsso, case_)
    case_ == "sweptwing_flat" ||
        error("mesh: only 'sweptwing_flat' is wired for now (got '$case_').")
    src = joinpath(ROOT, "PreProcessing", "modules", "directFlatPlateModule",
                   "postProcessing", "midPlane.csv")
    isfile(src) ||
        error("mesh: source midPlane not found:\n  $src\n(run PreProcessing DirectFlatPlate first)")
    dst = joinpath(dehnsso, "baseflow", "output", "benchmark", "bf_$(case_).csv")
    @info "airPower ▶ instAbility DeHNSSo mesh" src dst
    cp(src, dst; force = true)
    gridgen = joinpath(dehnsso, "gridgen", "benchmark", "$(case_).m")
    isfile(gridgen) || error("mesh: gridgen script not found: $gridgen")
    matlab_batch("run('$gridgen')")
end

# run: run the solver caller, then save StabRes+StabGrid for PostProcessing.
# The filename encodes the run's major inputs — <case>_<lin|nl>_N<N>_lz<λz>mm_A<A0> —
# built MATLAB-side (Stab/Opt live in the caller, not inputs.jl). In the name,
# '.'->'d' and '-'->'m' (e.g. 7.5->7d5, 5e-05->5em05) to avoid dots/dashes.
function dehnsso_run(dehnsso, case_)
    caller = joinpath(dehnsso, "callers", "benchmark", "$(case_).m")
    isfile(caller) || error("run: caller not found: $caller")
    fieldsdir = joinpath(ROOT, "PostProcessing", "io", "input")
    mkpath(fieldsdir)
    @info "airPower ▶ instAbility DeHNSSo run" caller fieldsdir
    stmt = string(
        "run('$caller'); ",
        "lz = 2*pi*StabGrid.lref/Stab.beta_0*1e3; ",
        "lin = 'nl'; if strcmpi(Opt.linear,'on'), lin = 'lin'; end; ",
        "nm = sprintf('$(case_)_%s_N%d_lz%.4gmm_A%.3g', lin, Stab.N, lz, Stab.A0_fund); ",
        "nm = strrep(strrep(nm, '.', 'd'), '-', 'm'); ",
        "out = fullfile('$fieldsdir', [nm '.mat']); ",
        "save(out, 'StabRes', 'StabGrid'); ",
        "fprintf('\\n[airPower] saved -> %s\\n', out)")
    matlab_batch(stmt)
end

# write the PostProcessing block of inputs.jl to the generated MATLAB config
function write_pp_config(pp, task)
    gen = joinpath(ROOT, "PostProcessing", "inputs_gen.m")
    open(gen, "w") do io
        println(io, "% AUTO-GENERATED by airPower dispatcher (run.jl) from inputs.jl — do not edit")
        println(io, "inp.airPowerRoot   = $(mat_lit(ROOT));")
        println(io, "inp.task           = $(mat_lit(task));")
        println(io, "inp.loadMode       = $(mat_lit(pp.loadMode));")
        println(io, "inp.caseType       = $(mat_lit(pp.caseType));")
        println(io, "inp.fieldsFile     = $(mat_lit(pp.fieldsFile));")
        println(io, "inp.modeIdx        = $(mat_lit(pp.modeIdx));")
        println(io, "inp.plotUProfiles  = $(mat_lit(pp.plotUProfiles));")
        println(io, "inp.validation     = $(mat_lit(pp.validation));")
        println(io, "inp.valXcZoom      = $(mat_lit(pp.valXcZoom));")
        println(io, "inp.valShareX      = $(mat_lit(pp.valShareX));")
        println(io, "inp.valYTop        = $(mat_lit(pp.valYTop));")
        # pulled from the PreProcessing blocks (no PostProcessing duplication):
        println(io, "inp.valPIV         = $(mat_lit(inp.VAL.valPIV));")
        println(io, "inp.valGen         = $(mat_lit(inp.VAL.Gen));")
        println(io, "inp.valCase        = $(mat_lit(inp.VAL.Case));")
        println(io, "inp.xInlet         = $(mat_lit(inp.DFP.xInlet));")
        # TTCP airfoil geometry (from the TunnelCase block) — needed by the
        # PostProcessing validation to map PIV x/c stations onto the curved-wall
        # stability grid (inverse rigid transform: physical wall point -> x/c).
        # Written for every case; only consumed by the TTCP validation branch.
        println(io, "inp.airfoilChord    = $(mat_lit(inp.TTCP.tunnel.chord));")     # [m]
        println(io, "inp.airfoilAlphaDeg = $(mat_lit(inp.TTCP.tunnel.alphaDeg));")  # [deg]
        println(io, "inp.airfoilXCenter  = $(mat_lit(inp.TTCP.tunnel.xCenter));")   # [m]
        println(io, "inp.airfoilYCenter  = $(mat_lit(inp.TTCP.tunnel.yCenter));")   # [m]
        println(io, "inp.ro.loadAnalysis = $(mat_lit(pp.ro.loadAnalysis));")
        # Shared plot window. An inputs.jl predating the `plot` block still runs:
        # bufferFrac/yWallFrac carry over from `ro`. The y-ceilings do not — the old
        # single `ro.yMax` could not say which load mode it meant, so it is left to
        # the yWallFrac fallback rather than guessed at.
        plt = hasproperty(pp, :plot) ? pp.plot : (;)
        if !hasproperty(pp, :plot)
            @warn "inputs.jl has no PostProcessing.plot block — plot windows fall back " *
                  "to yWallFrac (30% of the domain). Move bufferFrac/yMax out of `ro` " *
                  "into `plot` as yMaxFields [delta_0] and yMaxBF [m]."
        end
        getk(nt, k, d) = hasproperty(nt, k) ? getproperty(nt, k) :
                         (hasproperty(pp.ro, k) ? getproperty(pp.ro, k) : d)
        println(io, "inp.plot.bufferFrac = $(mat_lit(getk(plt, :bufferFrac, 0.85)));")
        println(io, "inp.plot.yWallFrac  = $(mat_lit(getk(plt, :yWallFrac, 0.30)));")
        println(io, "inp.plot.yMaxFields = $(mat_lit(hasproperty(plt, :yMaxFields) ? plt.yMaxFields : []));")
        println(io, "inp.plot.yMaxBF     = $(mat_lit(hasproperty(plt, :yMaxBF)     ? plt.yMaxBF     : []));")
    end
    return gen
end

# --- PostProcessing: export inp.PostProcessing -> generated MATLAB config, run.m ---
# Special task `config` regenerates the MATLAB config from inputs.jl WITHOUT running
# MATLAB — used by run.m to refresh a stale config in the direct-MATLAB workflow.
function run_postprocessing(args)
    isempty(args) && error("PostProcessing needs a task: " *
                           "julia run.jl PostProcessing <importData|reynoldsOrrProdTerms|config>")
    pp  = inp.PostProcessing
    sub = args[1]
    if sub == "config"
        gen = write_pp_config(pp, pp.task)
        @info "airPower ▶ PostProcessing config regenerated from inputs.jl" task = pp.task gen
        return
    end
    gen = write_pp_config(pp, sub)
    @info "airPower ▶ PostProcessing (MATLAB)" task = sub gen
    # AIRPOWER_PP_DISPATCH signals run.m that this config is authoritative (carries
    # the CLI task) — so run.m loads it as-is instead of re-syncing from inputs.jl.
    stmt = "run('$(joinpath(ROOT, "PostProcessing", "run.m"))')"
    matlab_batch(stmt; AIRPOWER_PP_DISPATCH = "1")
end

# --- main ---
if isempty(ARGS)
    usage()
    exit(0)
end

stage = ARGS[1]
rest  = ARGS[2:end]

if stage == "PreProcessing"
    run_preprocessing(rest)
elseif stage == "instAbility"
    run_instability(rest)
elseif stage == "PostProcessing"
    run_postprocessing(rest)
else
    @error "Unknown stage" stage
    usage()
    exit(1)
end
