using DelimitedFiles, LaTeXStrings, MAT

# Shared styling for the five U_e (freestream-edge) methods.
const METHOD_COLORS = (vortT=:royalblue, vortM=:darkorange, max=:firebrick,
                       fix=:forestgreen, bern=:purple)
const METHOD_LABELS = (vortT="vorticityIntegralTrapezoidal", vortM="vorticityIntegralMidpoint",
                       max="maxProfile", fix="fixedHeight", bern="pressureBernoulli")

# Caveat shown on the U_e comparison figures: pressureBernoulli returns the
# in-plane SPEED √(u²+v²) (total-head recovery), whereas the other four return
# the tangential COMPONENT u_t = u·t̂. They differ by the wall-normal velocity,
# most noticeably near the high-curvature leading edge. (Two short lines — GR
# mathtext has no \text, so wording avoids hyphenated words in math mode.)
const UE_DEF_NOTE1 = L"\mathrm{pressureBernoulli:}\ U_e=\sqrt{u^2+v^2}\ \ (\mathrm{in\ plane\ speed})"
const UE_DEF_NOTE2 = L"\mathrm{other\ methods:}\ U_e=u_t\ \ (\mathrm{tangential\ component})"

"""
    plot_bl_metrics(case_path; savedir)

Read `postProcessing/BLQuantities.csv` (single configured method) and
produce a 2×3 figure of U_e, δ99, δ* versus the surface arclength S (top
row) and versus the chord fraction xi = x/c (bottom row):

    ┌──────────┬──────────┬──────────┐
    │  U_e(S)  │ δ99(S)   │  δ*(S)   │
    ├──────────┼──────────┼──────────┤
    │ U_e(xi)  │ δ99(xi)  │  δ*(xi)  │
    └──────────┴──────────┴──────────┘

Saves to `<savedir>/BLQuantities<caseLabel>.png` and returns the figure.
"""
function plot_bl_metrics(case_path::AbstractString;
                         savedir::AbstractString=case_path)
    csv_path = joinpath(case_path, "postProcessing", "BLQuantities.csv")
    if !isfile(csv_path)
        @warn "BLQuantities.csv not found at $csv_path"
        return nothing
    end

    raw = readdlm(csv_path, ','; skipstart=1)
    # Columns: x, s, xi, Ue, d99, dstar, Theta
    S   = Float64.(raw[:, 2])
    xi  = Float64.(raw[:, 3])
    Ue  = Float64.(raw[:, 4])
    d99 = Float64.(raw[:, 5])
    dst = Float64.(raw[:, 6])

    @info "BL metrics (production): $(length(S)) wall faces"
    mkpath(savedir)

    common = (
        framestyle    = :box,
        grid          = true,
        gridalpha     = 0.3,
        tickfontsize  = 10,
        guidefontsize = 12,
        titlefontsize = 13,
        left_margin   = 8Plots.mm,
        bottom_margin = 6Plots.mm,
        right_margin  = 4Plots.mm,
        top_margin    = 4Plots.mm,
        linewidth     = 2,
        dpi           = 200,
        legend        = false,
    )

    function panel(x, y; xlabel, ylabel, title, yscale=1.0, color)
        plot(x, yscale .* y;
             xlabel=xlabel, ylabel=ylabel, title=title, color=color, common...)
    end

    pUe_S   = panel(S,  Ue;             xlabel=L"S\ \mathrm{[m]}",
                    ylabel=L"U_e\ \mathrm{[m/s]}",          title="Edge velocity",
                    color=:royalblue)
    pd99_S  = panel(S,  d99; yscale=1e3, xlabel=L"S\ \mathrm{[m]}",
                    ylabel=L"\delta_{99}\ \mathrm{[mm]}",   title=L"\delta_{99}",
                    color=:royalblue)
    pdst_S  = panel(S,  dst; yscale=1e3, xlabel=L"S\ \mathrm{[m]}",
                    ylabel=L"\delta^*\ \mathrm{[mm]}",      title=L"\delta^*",
                    color=:royalblue)
    pUe_xi  = panel(xi, Ue;             xlabel=L"x/c",
                    ylabel=L"U_e\ \mathrm{[m/s]}",          title="",
                    color=:firebrick)
    pd99_xi = panel(xi, d99; yscale=1e3, xlabel=L"x/c",
                    ylabel=L"\delta_{99}\ \mathrm{[mm]}",   title="",
                    color=:firebrick)
    pdst_xi = panel(xi, dst; yscale=1e3, xlabel=L"x/c",
                    ylabel=L"\delta^*\ \mathrm{[mm]}",      title="",
                    color=:firebrick)

    fig = plot(pUe_S, pd99_S, pdst_S, pUe_xi, pd99_xi, pdst_xi;
               layout=(2, 3), size=(1500, 800))

    label = basename(case_path)
    outfile = joinpath(savedir, "BLQuantities$(label).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end

"""
    plot_bl_metrics_comparison(case_path; savedir, ref_mat)

Produce a 2×2 figure overlaying the five U_e methods (vorticityIntegralTrapezoidal,
vorticityIntegralMidpoint, maxProfile, fixedHeight, pressureBernoulli) on each
panel:

    ┌─────────────────┬─────────────────┐
    │     U_e(s)      │    δ99(s)       │
    ├─────────────────┼─────────────────┤
    │     δ*(s)       │     θ(s)        │
    └─────────────────┴─────────────────┘

The per-method metrics are computed in-memory via `compute_bl_all_methods`
(no intermediate `BLQuantities_compare.csv` is written). Saves to
`<savedir>/BLQuantitiesCompare<caseLabel>.png` and returns the figure.

For a direct U_e-vs-x/c overlay against the external `.mat` reference (no
arclength mapping), see `plot_ue_methods_xc`.
"""
function plot_bl_metrics_comparison(case_path::AbstractString;
                                    savedir::AbstractString=case_path)
    r = try
        compute_bl_all_methods(case_path)
    catch err
        @warn "Freestream-method comparison skipped — could not compute BL metrics" exception=err
        return nothing
    end

    # Sort wall faces by xi for a clean monotone line, then bundle per method.
    ord   = r.perm
    S     = r.ss[ord]
    Ue    = (vortT=r.Ue_vT[ord],  vortM=r.Ue_vM[ord],  max=r.Ue_m[ord],  fix=r.Ue_f[ord],  bern=r.Ue_b[ord])
    d99   = (vortT=r.d99_vT[ord], vortM=r.d99_vM[ord], max=r.d99_m[ord], fix=r.d99_f[ord], bern=r.d99_b[ord])
    dst   = (vortT=r.dst_vT[ord], vortM=r.dst_vM[ord], max=r.dst_m[ord], fix=r.dst_f[ord], bern=r.dst_b[ord])
    Theta = (vortT=r.th_vT[ord],  vortM=r.th_vM[ord],  max=r.th_m[ord],  fix=r.th_f[ord],  bern=r.th_b[ord])

    @info "BL metrics comparison: $(length(S)) wall faces × 5 methods"

    mkpath(savedir)

    common = (
        framestyle    = :box,
        grid          = true,
        gridalpha     = 0.3,
        tickfontsize  = 10,
        guidefontsize = 12,
        titlefontsize = 13,
        left_margin   = 8Plots.mm,
        bottom_margin = 6Plots.mm,
        right_margin  = 4Plots.mm,
        top_margin    = 4Plots.mm,
        linewidth     = 2,
        dpi           = 200,
    )

    colors = METHOD_COLORS
    labels = METHOD_LABELS

    function fivepanel(yv; ylabel, title, yscale=1.0, ylog=false)
        p = plot(S, yscale .* yv.vortT; label=labels.vortT, color=colors.vortT,
                 xlabel=L"S\ \mathrm{[m]}", ylabel=ylabel, title=title,
                 yscale = ylog ? :log10 : :identity, common...)
        plot!(p, S, yscale .* yv.vortM; label=labels.vortM, color=colors.vortM, common...)
        plot!(p, S, yscale .* yv.max;   label=labels.max,   color=colors.max,   common...)
        plot!(p, S, yscale .* yv.fix;   label=labels.fix,   color=colors.fix,   common...)
        plot!(p, S, yscale .* yv.bern;  label=labels.bern,  color=colors.bern,  common...)
        return p
    end

    p1 = fivepanel(Ue;    ylabel=L"U_e\ \mathrm{[m/s]}",     title="Edge velocity")
    # U_e-definition caveat. Bottom-right is the empty region here (the legend
    # occupies the top-left, and U_e rises with S so the lower-right is clear).
    xR  = maximum(S) - 0.01*(maximum(S) - minimum(S))
    yh  = maximum(filter(isfinite, Ue.vortM)); yl = minimum(filter(isfinite, Ue.bern))
    annotate!(p1, xR, yl + 0.12*(yh - yl), text(UE_DEF_NOTE1, 7, :gray35, :right, :bottom))
    annotate!(p1, xR, yl + 0.05*(yh - yl), text(UE_DEF_NOTE2, 7, :gray35, :right, :bottom))
    p2 = fivepanel(d99;   yscale=1e3, ylabel=L"\delta_{99}\ \mathrm{[mm]}", title=L"\delta_{99}")
    p3 = fivepanel(dst;   yscale=1e3, ylabel=L"\delta^*\ \mathrm{[mm]}",    title=L"\delta^*")
    p4 = fivepanel(Theta; yscale=1e3, ylabel=L"\theta\ \mathrm{[mm]}",      title=L"\theta")

    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1400, 900))

    label = basename(case_path)
    outfile = joinpath(savedir, "BLQuantitiesCompare$(label).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end

"""
    plot_ue_methods_xc(case_path; savedir, ref_mat)

Single-panel edge velocity U_e versus x/c, overlaying all five CFD freestream
methods (computed in-memory via `compute_bl_all_methods`) and, when `ref_mat` is
given, the external reference U_e from the flow-data `.mat` (`BL["x"]/BL["c"]`,
`BL["Ue"]`). Both the CFD and the reference are native in x/c, so the comparison
is direct — no arclength mapping is involved. Saves to
`<savedir>/UeMethodsXC<caseLabel>.png` and returns the figure.
"""
function plot_ue_methods_xc(case_path::AbstractString;
                            savedir::AbstractString=case_path,
                            ref_mat::Union{AbstractString, Nothing}=nothing)
    r = try
        compute_bl_all_methods(case_path)
    catch err
        @warn "Ue-vs-x/c comparison skipped — could not compute BL metrics" exception=err
        return nothing
    end
    ord = r.perm
    xi  = r.xis[ord]
    mkpath(savedir)

    common = (framestyle=:box, grid=true, gridalpha=0.3, linewidth=2,
              tickfontsize=10, guidefontsize=12, titlefontsize=13,
              left_margin=8Plots.mm, bottom_margin=6Plots.mm, dpi=200)

    p = plot(xi, r.Ue_vT[ord]; label=METHOD_LABELS.vortT, color=METHOD_COLORS.vortT,
             xlabel=L"x/c", ylabel=L"U_e\ \mathrm{[m/s]}", title="Edge velocity",
             legend=:bottomright, size=(900, 560), common...)
    plot!(p, xi, r.Ue_vM[ord]; label=METHOD_LABELS.vortM, color=METHOD_COLORS.vortM, common...)
    plot!(p, xi, r.Ue_m[ord];  label=METHOD_LABELS.max,   color=METHOD_COLORS.max,   common...)
    plot!(p, xi, r.Ue_f[ord];  label=METHOD_LABELS.fix,   color=METHOD_COLORS.fix,   common...)
    plot!(p, xi, r.Ue_b[ord];  label=METHOD_LABELS.bern,  color=METHOD_COLORS.bern,  common...)

    if ref_mat !== nothing
        if !isfile(ref_mat)
            @warn "valUe is on but the reference flow-data file was not found — \
                   reference Ue overlay skipped (check inp.VAL.externalToScaling.\
                   flowDataFile, incl. the .mat extension)" ref_mat
        else
            BL   = matread(ref_mat)["BL"]
            cref = BL["c"]; cref = cref isa Number ? Float64(cref) : Float64(first(cref))
            xr   = vec(Float64.(BL["x"])) ./ cref
            ur   = vec(Float64.(BL["Ue"]))
            o    = sortperm(xr); xr = xr[o]; ur = ur[o]
            keep = (xr .>= minimum(xi)) .& (xr .<= maximum(xi))   # clip to CFD x/c range
            any(keep) && plot!(p, xr[keep], ur[keep];
                  label="Reference (.mat)", color=:black, linestyle=:dash,
                  linewidth=2, marker=:circle, markersize=3, markerstrokewidth=0)
        end
    end

    # Flag the U_e definition difference (Bernoulli = in-plane speed; others =
    # tangential). Two lines anchored top-left, where U_e(x/c) is lowest → no
    # data overlap.
    x0  = minimum(xi) + 0.01*(maximum(xi) - minimum(xi))
    yhi = maximum(filter(isfinite, r.Ue_vM[ord]))
    ylo = minimum(filter(isfinite, r.Ue_b[ord]))
    dy  = yhi - ylo
    annotate!(p, x0, yhi,             text(UE_DEF_NOTE1, 7, :gray35, :left, :top))
    annotate!(p, x0, yhi - 0.06*dy,   text(UE_DEF_NOTE2, 7, :gray35, :left, :top))

    outfile = joinpath(savedir, "UeMethodsXC$(basename(case_path)).png")
    savefig(p, outfile)
    @info "Saved: $outfile"
    return p
end
