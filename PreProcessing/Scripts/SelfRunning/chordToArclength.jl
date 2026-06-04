#!/usr/bin/env julia
#
# Convert chord percentage to arclength S, with two modes:
#
#   DFP   - read .mat from InputOutput/AirfoilFlowData/
#           (file picked by inp.VAL.externalToScaling.flowDataFile).
#           Uses the precomputed BL.x, BL.S, BL.c table — same data the
#           DFP solver consumes, so the resulting S is consistent with
#           what DFP sees internally.
#           Output adds xCenter = S - inp.DFP.xInlet (DFP-domain coord).
#
#   TTCP  - read .dat from InputOutput/AirfoilGeometryData/
#           (file picked by inp.TTCP.tunnel.airfoilFile). Integrates the
#           upper-surface arclength on the fly from the Selig-format
#           airfoil coordinates, scaled by inp.TTCP.tunnel.chord. This
#           matches the arclength TTCP/AirfoilLECase's controlDict
#           computes internally.
#
# Usage:
#   julia chordToArclength.jl DFP 15            # one chord% (DFP mode)
#   julia chordToArclength.jl TTCP 10 15 20     # multiple (TTCP mode)
#
# .dat format note: Selig only (one (x,y) pair per line; starts at TE
# upper, walks counterclockwise around the LE back to TE lower). All
# airfoils in this project (M3J, NACA*) are in this format.

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

using Printf

function usage()
    println("Usage: julia chordToArclength.jl <DFP|TTCP> <chord%> [chord%] ...")
    println("  e.g. julia chordToArclength.jl DFP 15")
    println("       julia chordToArclength.jl TTCP 10 20 30")
end

if length(ARGS) < 2 || !(uppercase(ARGS[1]) in ("DFP", "TTCP"))
    usage(); exit()
end

mode        = uppercase(ARGS[1])
percentages = parse.(Float64, ARGS[2:end])

# --- helper: 1-D linear interpolation over a monotonically increasing table ---
function interp_S(x_ref, S_ref, xq)
    idx = clamp(searchsortedlast(x_ref, xq), 1, length(x_ref) - 1)
    t   = (xq - x_ref[idx]) / (x_ref[idx+1] - x_ref[idx])
    return (1 - t) * S_ref[idx] + t * S_ref[idx+1]
end

if mode == "DFP"
    using MAT
    flow_dir = joinpath(ROOT, "PreProcessing", "InputOutput", "AirfoilFlowData")
    fname    = inp.VAL.externalToScaling.flowDataFile
    mat_path = joinpath(flow_dir, fname)
    isfile(mat_path) || error(
        ".mat file not found: $mat_path\n" *
        "(set by inp.VAL.externalToScaling.flowDataFile)")

    BL     = matread(mat_path)["BL"]
    x_ref  = vec(BL["x"])
    S_ref  = vec(BL["S"])
    c_ref  = BL["c"]
    xInlet = inp.DFP.xInlet

    println("─────────────────────────────────────────────────────")
    println("  Chord [%]    x [m]        S [m]        xCenter [m]")
    println("                                         (S - xInlet)")
    println("─────────────────────────────────────────────────────")
    for pct in percentages
        x_chord = pct / 100.0 * c_ref
        S_val   = interp_S(x_ref, S_ref, x_chord)
        xCenter = S_val - xInlet
        @printf("  %6.1f      %8.5f     %8.5f     %8.5f\n",
                pct, x_chord, S_val, xCenter)
    end
    println("─────────────────────────────────────────────────────")
    @printf("Source: %s  (chord = %.5g m)\n", fname, c_ref)
    println("Note: xCenter is in DFP domain coordinates (x=0 at inlet).")

else  # TTCP
    geom_dir = joinpath(ROOT, "PreProcessing", "InputOutput", "AirfoilGeometryData")
    fname    = inp.TTCP.tunnel.airfoilFile
    dat_path = joinpath(geom_dir, fname)
    isfile(dat_path) || error(
        ".dat file not found: $dat_path\n" *
        "(set by inp.TTCP.tunnel.airfoilFile)")

    # Parse (x, y) pairs from the .dat; skip non-numeric lines and
    # any leading header (e.g. Lednicer's "n_upper n_lower" gets
    # dropped because its values are outside the normalised range).
    coords = NTuple{2,Float64}[]
    for line in readlines(dat_path)
        parts = split(strip(line))
        length(parts) >= 2 || continue
        a = tryparse(Float64, parts[1])
        b = tryparse(Float64, parts[2])
        (a === nothing || b === nothing) && continue
        # Drop obvious header lines (airfoil coords are normalised to [0,1])
        (a < -0.1 || a > 1.1 || abs(b) > 0.5) && continue
        push!(coords, (a, b))
    end
    isempty(coords) && error("No coordinate pairs parsed from $dat_path")

    # Selig format: starts at TE upper, walks counterclockwise around
    # LE to TE lower. Upper surface = points before LE, reversed so
    # they run LE → TE.
    le_idx = argmin(getindex.(coords, 1))
    upper  = coords[le_idx:-1:1]
    n_pts  = length(upper)
    n_pts < 3 && error("Upper surface has only $n_pts points — not Selig format?")

    chord = inp.TTCP.tunnel.chord
    xs_n  = getindex.(upper, 1)            # normalised xi ∈ [0,1]
    ys_n  = getindex.(upper, 2)
    S_n   = zeros(n_pts)                   # normalised arclength
    for k in 2:n_pts
        S_n[k] = S_n[k-1] + sqrt((xs_n[k]-xs_n[k-1])^2 + (ys_n[k]-ys_n[k-1])^2)
    end
    x_ref = xs_n .* chord                  # physical x [m]
    S_ref = S_n  .* chord                  # physical S [m]

    println("─────────────────────────────────────")
    println("  Chord [%]    x [m]        S [m]")
    println("─────────────────────────────────────")
    for pct in percentages
        x_chord = pct / 100.0 * chord
        S_val   = interp_S(x_ref, S_ref, x_chord)
        @printf("  %6.1f      %8.5f     %8.5f\n", pct, x_chord, S_val)
    end
    println("─────────────────────────────────────")
    @printf("Source: %s  (chord = %.5g m, n_upper = %d, S_max = %.5g m)\n",
            fname, chord, n_pts, S_ref[end])
end
