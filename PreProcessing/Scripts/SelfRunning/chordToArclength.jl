#!/usr/bin/env julia
#
# Convert chord percentage to arc-length S for bump placement.
#
# Usage:
#   julia chordToArclength.jl 15        # single value: 15% chord
#   julia chordToArclength.jl 10 15 20  # multiple values
#
# The arc-length S can be used as xCenter in inputs.jl DFP.wallModulation.
# For the DirectFlatPlate, xCenter = S - xInlet (domain coordinates).
#

# Find and load inputs
const ROOT = let
    dir = @__DIR__
    while !isfile(joinpath(dir, "inputs.jl"))
        parent = dirname(dir)
        parent == dir && error("inputs.jl not found above $(@__DIR__)")
        dir = parent
    end
    dir
end
include(joinpath(ROOT, "inputs.jl"))

using MAT, Printf

# Load reference data
flow_data_dir = joinpath(ROOT, "PreProcessing", "InputOutput", "AirfoilFlowData")
mat_files = filter(f -> endswith(f, ".mat"), readdir(flow_data_dir))
isempty(mat_files) && error("No .mat files found in $flow_data_dir")

mat_path = joinpath(flow_data_dir, mat_files[1])
BL = matread(mat_path)["BL"]
x_ref = vec(BL["x"])   # chordwise coordinate [m]
S_ref = vec(BL["S"])   # arc-length [m]
c_ref = BL["c"]        # chord [m]

xInlet = inp.DFP.xInlet

# Parse chord percentages from command line
if isempty(ARGS)
    println("Usage: julia chordToArclength.jl <chord%> [chord%] ...")
    println("Example: julia chordToArclength.jl 15")
    exit()
end

percentages = parse.(Float64, ARGS)

# Interpolate S at each chord percentage
println("─────────────────────────────────────────────────────")
println("  Chord [%]    x/c [m]      S [m]        xCenter [m]")
println("                                         (S - xInlet)")
println("─────────────────────────────────────────────────────")
for pct in percentages
    x_chord = pct / 100.0 * c_ref

    # Interpolate S at this x
    idx = searchsortedlast(x_ref, x_chord)
    idx = clamp(idx, 1, length(x_ref) - 1)
    t = (x_chord - x_ref[idx]) / (x_ref[idx+1] - x_ref[idx])
    S_val = (1 - t) * S_ref[idx] + t * S_ref[idx+1]

    xCenter = S_val - xInlet

    Printf.@printf("  %6.1f      %8.5f     %8.5f     %8.5f\n",
                    pct, x_chord, S_val, xCenter)
end
println("─────────────────────────────────────────────────────")
println("Note: xCenter is in domain coordinates (x=0 at inlet).")
