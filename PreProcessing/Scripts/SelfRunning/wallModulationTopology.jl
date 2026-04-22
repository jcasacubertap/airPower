#!/usr/bin/env julia
#
# Plot the wall modulation shape defined in inputs.jl
# Usage: julia wallModulationTopology.jl
#

# Find and load inputs.jl
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

using Plots
gr()
default(fontfamily = "Computer Modern")

# ── Smoothstep bump function ──────────────────────────────────────────
"""
    smoothstep(t, n)

Generalized smoothstep: S(0)=0, S(1)=1, S'(0)=S'(1)=0 for n≥2.
"""
smoothstep(t, n) = t^n / (t^n + (1 - t)^n)

"""
    wall_bump(x; A, xStart, xPeak, xEnd, p, q)

Evaluate the wall modulation height at position x.
Uses piecewise sigmoidal smoothstep for C2 continuity (p,q ≥ 3).
"""
function wall_bump(x; A, xStart, xPeak, xEnd, p, q)
    (x <= xStart || x >= xEnd) && return 0.0
    if x <= xPeak
        t = (x - xStart) / (xPeak - xStart)
        return A * smoothstep(t, p)
    else
        s = (x - xPeak) / (xEnd - xPeak)
        return A * (1.0 - smoothstep(s, q))
    end
end

# ── Read inputs and plot ──────────────────────────────────────────────
wm = inp.DFP.wallModulation

if !wm.enabled
    @info "wallModulation is disabled in inputs.jl — plotting with current parameters anyway"
end

if wm.mode != :single
    @warn "Only :single mode is supported for now"
    exit()
end

# Domain range
L = inp.DFP.domainLength
x = range(0, L, length=1000)
h = [wall_bump(xi; A=wm.A, xStart=wm.xStart, xPeak=wm.xPeak, xEnd=wm.xEnd, p=wm.p, q=wm.q) for xi in x]

fig = plot(x, h .* 1000;   # convert to mm for readability
    xlabel     = "x [m]",
    ylabel     = "h [mm]",
    title      = "Wall modulation (A=$(wm.A*1000) mm, p=$(wm.p), q=$(wm.q))",
    color      = :black,
    linewidth  = 2,
    legend     = false,
    framestyle = :box,
    grid       = true,
    gridalpha  = 0.3,
    tickfontsize   = 10,
    guidefontsize  = 12,
    titlefontsize  = 13,
    left_margin    = 8Plots.mm,
    bottom_margin  = 6Plots.mm,
    ylims = (-wm.A * 1000 * 0.5, wm.A * 1000 * 5),
    size       = (800, 400),
    dpi        = 200,
)

# Mark key positions
vline!(fig, [wm.xStart, wm.xPeak, wm.xEnd];
    color=:gray, linestyle=:dash, linewidth=0.8, label=false)

# Save
outdir = joinpath(ROOT, "PreProcessing", "InputOutput", "Plotting", "DirectFlatPlate")
mkpath(outdir)
outfile = joinpath(outdir, "wallModulationTopology.png")
savefig(fig, outfile)
@info "Saved: $outfile"
