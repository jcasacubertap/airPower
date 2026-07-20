#!/usr/bin/env julia
#
# airPower — PreProcessing orchestrator
#
# Usually invoked via the top-level airPower/run.jl dispatcher:
#   julia run.jl PreProcessing <module> <action>
# but also runnable directly:
#   julia PreProcessing/run.jl TunnelToCurvedPlate all        # full pipeline: tunnel → curved plate
#   julia PreProcessing/run.jl TunnelToCurvedPlate meshTunnel # mesh tunnel only
#   julia PreProcessing/run.jl TunnelToCurvedPlate runTunnel  # solve tunnel only
#   julia PreProcessing/run.jl TunnelToCurvedPlate map        # map tunnel → airfoil BCs
#   julia PreProcessing/run.jl TunnelToCurvedPlate runAirfoil # solve airfoil only
#   julia PreProcessing/run.jl DirectFlatPlate all            # direct flat-plate computation
#   julia PreProcessing/run.jl DirectFlatPlate viz            # visualization only
#

# This file lives in PreProcessing/, but ROOT is the airPower root so that every
# path below (inputs.jl, PreProcessing/src/..., module dirs) stays unchanged.
const ROOT = dirname(@__DIR__)

# Central inputs (single source of truth for all parameters)
include(joinpath(ROOT, "inputs.jl"))

# Include all source files
include(joinpath(ROOT, "PreProcessing", "src", "source", "backend.jl"))

# Auxiliary functions (must come before source files that use them)
include(joinpath(ROOT, "PreProcessing", "src", "auxiliary", "leastSquares.jl"))
include(joinpath(ROOT, "PreProcessing", "src", "auxiliary", "falknerSkan.jl"))

# Plotting & validation (must come before modules that reference plot functions)
using Plots, DelimitedFiles, Glob, Printf, LaTeXStrings, Statistics
default(fontfamily = "Computer Modern")
const PV_DIR = joinpath(ROOT, "PreProcessing", "src", "plottingAndValidation")
include(joinpath(PV_DIR, "residuals.jl"))
include(joinpath(PV_DIR, "fields.jl"))
include(joinpath(PV_DIR, "profiles.jl"))
include(joinpath(PV_DIR, "wallgeometry.jl"))
include(joinpath(PV_DIR, "wallquantities.jl"))
include(joinpath(PV_DIR, "blmetrics.jl"))
include(joinpath(PV_DIR, "experimentalvalidation.jl"))
include(joinpath(PV_DIR, "dfpvalidation.jl"))
include(joinpath(PV_DIR, "iblTTCPComparison.jl"))

# External-to-scaling preprocessing
include(joinpath(ROOT, "PreProcessing", "src", "source", "externalToScaling.jl"))

# BL-metrics computation (also run standalone in postAirfoil); included so the
# viz step can build the freestream-method comparison figure in-memory, with no
# intermediate BLQuantities_compare.csv written to disk.
include(joinpath(ROOT, "PreProcessing", "src", "source", "blMetrics.jl"))

# Third-party boundary-layer (IBL) solver bridge (glue to external/ solver)
include(joinpath(ROOT, "PreProcessing", "src", "thirdParty", "bridge", "runIBL.jl"))

# Wall-modulation scaling: resolve a bump's {A, xCenter} from {Re_k, A/δ*} on the
# IBL baseline (DFP). Included after runIBL (uses run_ibl_solver) and before the
# modules (backend/directFlatPlate use dfp_wm at run time).
include(joinpath(ROOT, "PreProcessing", "src", "source", "wallModulationScaling.jl"))

# Module definitions and pipeline
include(joinpath(ROOT, "PreProcessing", "src", "source", "pipeline.jl"))
include(joinpath(ROOT, "PreProcessing", "src", "source", "tunnelToCurvedPlate.jl"))
include(joinpath(ROOT, "PreProcessing", "src", "source", "directFlatPlate.jl"))

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
