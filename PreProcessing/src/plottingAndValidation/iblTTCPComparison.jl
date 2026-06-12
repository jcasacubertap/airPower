#
# iblTTCPComparison — δ99/δ* validation of the spectral IBL solver against the
# TTCP airfoilLE CFD, plotted vs x/c (axis B).
#
# The IBL is run with the TTCP chosen-method edge velocity (its BLQuantities Ue),
# and its inflow virtual origin is δ*-calibrated to the TTCP inlet (option B), so
# both share an identical edge velocity and the δ* curves start coincident. The
# δ99/δ* integral procedure is `compute_bl_integrals` from blMetrics.jl — the
# SAME one used for the TTCP CFD — applied to the IBL with its (invariant) top-
# row freestream as Ue.
#
# This is .mat-free: it uses only the TTCP simulation output (BLQuantities.csv).

if !@isdefined(compute_bl_integrals)
    include(joinpath(@__DIR__, "..", "source", "blMetrics.jl"))
end

using DelimitedFiles, LaTeXStrings

"Linear interpolation of `ys` (sampled at sorted `xs`) at `xq`, clamped at ends."
function _interp1(xs::AbstractVector, ys::AbstractVector, xq::Real)
    xq <= xs[1]   && return float(ys[1])
    xq >= xs[end] && return float(ys[end])
    i = searchsortedlast(xs, xq); i = clamp(i, 1, length(xs) - 1)
    t = (xq - xs[i]) / (xs[i+1] - xs[i])
    return (1 - t) * ys[i] + t * ys[i+1]
end

"""
    read_bl_quantities(path) → NamedTuple | nothing

Header-aware reader for a BLQuantities CSV. Returns a NamedTuple keyed by the
column names (e.g. `x, s, xi, Ue, d99, dstar, Theta`). `nothing` if absent.
"""
function read_bl_quantities(path::AbstractString)
    isfile(path) || return nothing
    data, hdr = readdlm(path, ',', Float64, '\n'; header=true)
    names = Tuple(Symbol.(strip.(vec(hdr))))
    return NamedTuple{names}(Tuple(data[:, j] for j in 1:length(names)))
end

"""
    compute_ibl_bl_metrics(bl; max_n) → NamedTuple

δ99/δ*/θ of the IBL solution per station. Ue is the top-row (freestream) value
of each Chebyshev profile — invariant in y, i.e. the prescribed edge velocity.
Returns `(x, Ue, d99, dstar, Theta)`; `x` is the streamwise offset from the IBL
inlet (`bl.Xgrid`).
"""
function compute_ibl_bl_metrics(bl; max_n::Real)
    n  = length(bl.Xgrid)
    x  = collect(Float64, bl.Xgrid)
    yv = collect(Float64, bl.Y)
    Ue = fill(NaN, n); d99 = fill(NaN, n); dstar = fill(NaN, n); th = fill(NaN, n)
    for j in 1:n
        uj    = Float64.(bl.u[:, j])
        Ue[j] = uj[1]                                   # top row = freestream
        d99[j], dstar[j], th[j] = compute_bl_integrals(yv, uj, Ue[j]; max_n=max_n)
    end
    return (x = x, Ue = Ue, d99 = d99, dstar = dstar, Theta = th)
end

"""
    run_ibl_ttcp(blq) → bl | nothing

Run the δ*-calibrated spectral IBL for the TTCP airfoilLE, driven by the TTCP
`BLQuantities` edge velocity. `blq` (from `read_bl_quantities`) must carry
`s, Ue, dstar`. The marching span is the arclength range `[s_inlet, s_outlet]`;
the inflow virtual origin is calibrated so the FSC δ* matches the inlet `dstar`.
"""
function run_ibl_ttcp(blq)
    a = inp.TTCP.airfoilLE; f = inp.TTCP.flow
    perm   = sortperm(blq.s)
    s      = blq.s[perm]; Ue_t = blq.Ue[perm]; dstar_t = blq.dstar[perm]
    s_in, s_out = s[1], s[end]
    Lspan  = s_out - s_in

    # Shared IBL grid (inp.VAL); H is the height of the exported data box.
    nx  = Int(inp.VAL.blnx)
    ny  = Int(inp.VAL.blNcheb)
    y_i = Float64(inp.VAL.blYi)
    H   = Float64(a.exportHeight)

    sg = range(s_in, s_out, length=nx)
    Ue = [_interp1(s, Ue_t, sq) for sq in sg]
    dstar_in = dstar_t[1]

    @info "TTCP IBL (δ*-calibrated)" nx=nx ny=ny H=H y_i=y_i s_in=round(s_in, digits=4) s_out=round(s_out, digits=4) dstar_in_mm=round(1e3*dstar_in, digits=4)

    return run_ibl(; S=0.0, Lspan=Lspan, H=H, nu=f.freeStreamViscosity,
                   We=f.freeStreamVelocitySpanwise, Ue=Ue, nx=nx, ny=ny, y_i=y_i,
                   calibrate_dstar=dstar_in)
end

"""
    plot_ttcp_ibl_comparison(case_path; savedir) → fig | nothing

:vizAirfoil step (valBL). Read the TTCP `BLQuantities.csv`, run the δ*-calibrated
IBL inline (no files left behind), compute its δ99/δ*, map stations to `x/c` via
the TTCP `s↔xi` relation, and plot δ99 and δ* vs `x/c` — TTCP (symbols) vs IBL
(dashed line). Everything is in memory; only the figure is saved.
"""
function plot_ttcp_ibl_comparison(case_path::AbstractString; savedir::AbstractString)
    ttcp = read_bl_quantities(joinpath(case_path, "postProcessing", "BLQuantities.csv"))
    if ttcp === nothing
        @warn "TTCP BLQuantities.csv not found — run :postAirfoil first; skipping IBL comparison."
        return nothing
    end

    bl = run_ibl_ttcp(ttcp)
    bl === nothing && return nothing
    m = compute_ibl_bl_metrics(bl; max_n=Float64(inp.TTCP.airfoilLE.exportHeight))

    # Map IBL streamwise offset x → physical arclength s → x/c (TTCP s↔xi).
    perm = sortperm(ttcp.s); s = ttcp.s[perm]; xi = ttcp.xi[perm]
    s_in = s[1]
    xi_ibl = [_interp1(s, xi, s_in + xx) for xx in m.x]
    ibl = (xi=xi_ibl, d99=m.d99, dstar=m.dstar)

    mkpath(savedir)
    common = (
        framestyle = :box, grid = true, gridalpha = 0.3,
        tickfontsize = 10, guidefontsize = 12, legendfontsize = 9,
        left_margin = 8Plots.mm, bottom_margin = 6Plots.mm,
        right_margin = 4Plots.mm, top_margin = 4Plots.mm, dpi = 200,
    )
    mk = (tval, ival, ylab) -> begin
        p = plot(; xlabel = L"x/c", ylabel = ylab, legend = :topleft, common...)
        plot!(p, ttcp.xi, 1e3 .* tval;
              label = "TTCP", color = :royalblue, seriestype = :scatter,
              markersize = 2.5, markerstrokewidth = 0)
        plot!(p, ibl.xi, 1e3 .* ival;
              label = "IBL", color = :firebrick, linestyle = :dash, linewidth = 2)
        return p
    end
    p1 = mk(ttcp.d99,   ibl.d99,   L"\delta_{99} \ \mathrm{[mm]}")
    p2 = mk(ttcp.dstar, ibl.dstar, L"\delta^{*} \ \mathrm{[mm]}")
    fig = plot(p1, p2; layout = (1, 2), size = (1150, 430))
    outfile = joinpath(savedir, "blMetricsComparison$(basename(case_path)).png")
    savefig(fig, outfile); @info "Saved: $outfile"
    return fig
end
