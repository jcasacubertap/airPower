#!/usr/bin/env julia
#
# aeroForces — integrate static pressure and wall shear stress along the
# (suction-side) wall over a user-specified chord-fraction range, and project
# the results onto the global x axis.
#
# Reads:
#   <case>/constant/inputParam                 (Upinlet + aeroForces sub-dicts)
#   <case>/postProcessing/wallQuantities.csv   (x, s, p, tau, yPlus)
#   <case>/postProcessing/wallCoordinates.csv  (x, y, used to recover xi)
#
# Writes:
#   <case>/postProcessing/aeroForces.csv
#     xiStart,xiEnd,F_p_avg,F_p_x,F_p_y,F_tau_x
#     <…>
#
# Quantities — all are FORCES ON THE WALL FROM THE FLUID
# (kinematic units, m³/s² ≡ N/m per span / ρ; multiply by ρ for N/m per span):
#
#   F_p_avg = ⟨p⟩·arclength = ∫ p dS
#       Scalar "mean-pressure load" on the slice — what the perpendicular
#       force would be on a flat wall of the same arclength under uniform
#       pressure equal to the slice average. Geometry-blind, depends only on
#       p and s. Sign:  F_p_avg > 0 → wall is pushed (compressed),
#                       F_p_avg < 0 → wall is sucked (ripped off).
#
#   F_p_x   = −∫ p·n_x dS                  pressure force on the wall, x-component
#   F_p_y   = −∫ p·n_y dS                  pressure force on the wall, y-component
#   F_tau_x = +∫ τ_w·t_x dS                friction force on the wall, x-component
#
# Sign convention for the projected force components:
#   • pressure force per area on the wall = −p·n̂   (n̂ outward; fluid pushes wall in)
#   • friction force per area on the wall = +τ_w·t̂ (t̂ tangent toward +xi)
#
# All face geometry helpers (read_subdict, build_geom, surface_frame,
# project_wall_to_xi) are reused from blMetrics.jl via include.
#
# Invocation:
#   julia aeroForces.jl <caseDir>    # caseDir defaults to "."
#
# Jordi Casacuberta, 2026

using DelimitedFiles
using Printf

# Reuse geometry + parser helpers from the BL post-processor
include(joinpath(@__DIR__, "blMetrics.jl"))

function main_aeroforces(case_dir::AbstractString)
    inputParam = joinpath(case_dir, "constant", "inputParam")
    wallQ      = joinpath(case_dir, "postProcessing", "wallQuantities.csv")
    wallCoords = joinpath(case_dir, "postProcessing", "wallCoordinates.csv")

    @printf "aeroForces: case = %s\n" abspath(case_dir)
    isfile(inputParam) || error("missing $(inputParam)")
    isfile(wallQ)      || error("missing $(wallQ) — run runPostProcess first")
    isfile(wallCoords) || error("missing $(wallCoords) — run runPostProcess first")

    up = read_subdict(inputParam, "Upinlet")
    af = read_subdict(inputParam, "aeroForces")
    xiS = af["xiStart"]
    xiE = af["xiEnd"]
    xiS < xiE || error("aeroForces: xiStart ($xiS) must be less than xiEnd ($xiE)")

    # Read wall quantities (x, s, p, tau, yPlus) and wall coords (x, y)
    wq, _ = readdlm(wallQ, ',', Float64, '\n'; header=true)
    x_w = wq[:, 1]; s_w = wq[:, 2]; p_w = wq[:, 3]; tau_w = wq[:, 4]
    nW = length(x_w)

    wc, _ = readdlm(wallCoords, ',', Float64, '\n'; header=true)
    y_w = wc[:, 2]
    length(y_w) == nW || error(
        "wallQuantities.csv ($(nW) rows) and wallCoordinates.csv ($(length(y_w)) rows) disagree")

    # Per-face chord fraction xi (via projection onto the refined surface)
    geom    = build_geom(up)
    xis, _  = project_wall_to_xi(geom, x_w, y_w)
    frames  = [surface_frame(geom, xis[j]) for j in 1:nW]

    # Filter to the requested chord-fraction range, sort by arclength s
    sel = findall(j -> (xiS <= xis[j] <= xiE), 1:nW)
    isempty(sel) && error(
        "aeroForces: no wall faces with xi ∈ [$xiS, $xiE]  (available xi range: $(minimum(xis)) … $(maximum(xis)))")
    sort!(sel; by = j -> s_w[j])
    nSel = length(sel)

    # Trapezoidal weights from the sorted s positions:
    # face k gets dS_k = ½·(s_{k+1} − s_{k−1}), end faces get half-distance.
    s_sel = s_w[sel]
    dS = Vector{Float64}(undef, nSel)
    @inbounds for k in 1:nSel
        s_lo = k == 1     ? s_sel[1]   : 0.5 * (s_sel[k-1] + s_sel[k])
        s_hi = k == nSel  ? s_sel[end] : 0.5 * (s_sel[k]   + s_sel[k+1])
        dS[k] = s_hi - s_lo
    end

    # Raw integrals (signed by geometry — not yet forces in the standard
    # sense; see sign convention in the file docstring).
    I_p   = 0.0                                  # ∫ p dS
    I_p_x = 0.0                                  # ∫ p·n_x dS
    I_p_y = 0.0                                  # ∫ p·n_y dS
    I_t_x = 0.0                                  # ∫ τ_w·t_x dS
    for k in 1:nSel
        j  = sel[k]
        nx, ny, tx, _ = frames[j]
        I_p   += p_w[j]           * dS[k]
        I_p_x += p_w[j]   * nx    * dS[k]
        I_p_y += p_w[j]   * ny    * dS[k]
        I_t_x += tau_w[j] * tx    * dS[k]
    end

    # Final force-like quantities on the wall:
    #   F_p_avg : mean-pressure load on the slice = ⟨p⟩·arclength = ∫p dS
    #             Diagnostic of net normal pressure on the section, geometry-
    #             blind. Sign: + = wall pushed,  − = wall sucked off.
    #   F_p_x   : x-component of the actual pressure force on the wall
    #             (per area, this is −p·n̂; we integrate and project on x).
    #   F_p_y   : y-component of the actual pressure force on the wall.
    #   F_tau_x : x-component of the actual friction force on the wall.
    F_p_avg = I_p
    F_p_x   = -I_p_x
    F_p_y   = -I_p_y
    F_tau_x = +I_t_x

    # Write the output file (with a short sign-convention preamble; readers
    # should be configured to skip lines beginning with '#').
    outFile = joinpath(case_dir, "postProcessing", "aeroForces.csv")
    open(outFile, "w") do io
        println(io, "# aeroForces — forces on the wall from the fluid, integrated over")
        println(io, "# xi ∈ [xiStart, xiEnd] along the upper-surface arclength.")
        println(io, "# Kinematic units (m^3/s^2); multiply by ρ for N/m per span.")
        println(io, "# Sign convention:")
        println(io, "#   F_p_avg = ⟨p⟩·arclength = ∫ p dS")
        println(io, "#             > 0 → wall pushed (net positive pressure)")
        println(io, "#             < 0 → wall sucked off (net suction)")
        println(io, "#   F_p_x   = −∫ p·n_x dS   pressure force on the wall, x-component")
        println(io, "#   F_p_y   = −∫ p·n_y dS   pressure force on the wall, y-component")
        println(io, "#                            (n̂ is the outward wall normal; fluid")
        println(io, "#                             pushes the wall in the −n̂ direction)")
        println(io, "#   F_tau_x = +∫ τ_w·t_x dS friction force on the wall, x-component")
        println(io, "#                            (t̂ is the wall tangent toward +xi)")
        println(io, "xiStart,xiEnd,F_p_avg,F_p_x,F_p_y,F_tau_x")
        @printf(io, "%.6g,%.6g,%.10g,%.10g,%.10g,%.10g\n",
                xiS, xiE, F_p_avg, F_p_x, F_p_y, F_tau_x)
    end

    # Summary
    arc_total = sum(dS)
    @printf "  xi window = [%.4f, %.4f]   (%d wall faces)\n" xiS xiE nSel
    @printf "  s window  = [%.6f, %.6f] m   total arclength = %.6f m\n" s_sel[1] s_sel[end] arc_total
    @printf "  Forces on the wall (kinematic, multiply by ρ for N/m per span):\n"
    @printf "    F_p_avg = %+.6g  [m^3/s^2]   (= ⟨p⟩·arclength = %+.6g · %.6f)\n" F_p_avg (F_p_avg/arc_total) arc_total
    @printf "    F_p_x   = %+.6g              (pressure force, x-component)\n" F_p_x
    @printf "    F_p_y   = %+.6g              (pressure force, y-component)\n" F_p_y
    @printf "    F_tau_x = %+.6g              (friction force, x-component)\n" F_tau_x
    @printf "  → %s\n" outFile

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    case_dir = length(ARGS) >= 1 ? ARGS[1] : "."
    main_aeroforces(case_dir)
end
