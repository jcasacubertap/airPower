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
include(joinpath(ROOT, "PreProcessing", "Scripts", "Source", "backend.jl"))

using Plots
gr()
default(fontfamily = "Computer Modern")

# ── Read inputs ───────────────────────────────────────────────────────
wm = merge(inp.wallModulation, inp.DFP.wallModulation)

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
h = [wall_bump(xi, wm) for xi in x]

# Determine bump extent and title for either shape
if wm.shape == :sigmoidal
    bump_xs = wm.xStart
    bump_xe = wm.xEnd
    bump_xp = wm.xPeak
    title_str = "Sigmoidal (A=$(wm.A*1000) mm, p=$(wm.p), q=$(wm.q))"
elseif wm.shape == :esn
    bump_xs, bump_xe = _esn_geometry(wm)
    bump_xp = wm.xCenter
    title_str = "ESN (A=$(wm.A*1000) mm, ε=$(wm.epsilon), R=$(wm.R))"
end

fig = plot(x, h .* 1000;   # convert to mm
    xlabel     = "x [m]",
    ylabel     = "h [mm]",
    title      = title_str,
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
vline!(fig, [bump_xs, bump_xp, bump_xe];
    color=:gray, linestyle=:dash, linewidth=0.8, label=false)

# Save
outdir = joinpath(ROOT, "PreProcessing", "InputOutput", "Plotting", "DirectFlatPlate")
mkpath(outdir)
outfile = joinpath(outdir, "wallModulationTopology.png")
savefig(fig, outfile)
@info "Saved: $outfile"
