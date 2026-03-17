#!/usr/bin/env julia
#
# airPower orchestrator
#
# Usage:
#   julia run.jl TunnelToCurvedPlate all   # two-step: tunnel → curved plate
#   julia run.jl DirectFlatPlate all       # direct flat-plate computation
#   julia run.jl DirectFlatPlate viz       # visualization only
#

const ROOT = @__DIR__

# Central inputs (single source of truth for all parameters)
include(joinpath(ROOT, "inputs.jl"))

# Include all source files
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "backend.jl"))

# Plotting (must come before modules that reference plot functions)
using Plots, DelimitedFiles, Glob, Printf, LaTeXStrings
default(fontfamily = "Computer Modern")
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "residuals.jl"))
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "fields.jl"))
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "profiles.jl"))

# Module definitions and pipeline
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "pipeline.jl"))
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "TunnelToCurvedPlate.jl"))
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "DirectFlatPlate.jl"))

# --- Main ---
# Only run automatically when called from the command line (julia run.jl ...),
# not when loaded interactively via include("run.jl") in the REPL.
if !isempty(ARGS)
    module_name = ARGS[1]
    action      = length(ARGS) >= 2 ? ARGS[2] : "all"

    @info "airPower orchestrator" mod=module_name action=action
    run_module(module_name, action; root=ROOT)
else
    @info "airPower loaded. Usage: run_module(\"TunnelToCurvedPlate\", \"viz\"; root=ROOT)"
end
