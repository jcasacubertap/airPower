using DelimitedFiles, LaTeXStrings

"""
    plot_bl_metrics(case_path; savedir)

Read `postProcessing/BLQuantities_compare.csv` and produce a 2×2 figure
overlaying the three U_e methods (vorticityIntegral, maxProfile, fixedHeight)
on each panel:

    ┌─────────────────┬─────────────────┐
    │     U_e(s)      │    δ99(s)       │
    ├─────────────────┼─────────────────┤
    │     δ*(s)       │     θ(s)        │
    └─────────────────┴─────────────────┘

Saves to `<savedir>/BLQuantities<caseLabel>.png` and returns the figure.
"""
function plot_bl_metrics(case_path::AbstractString;
                         savedir::AbstractString=case_path)
    csv_path = joinpath(case_path, "postProcessing", "BLQuantities_compare.csv")
    if !isfile(csv_path)
        @warn "BLQuantities_compare.csv not found at $csv_path"
        return nothing
    end

    raw = readdlm(csv_path, ','; skipstart=1)
    S   = Float64.(raw[:, 2])
    # Columns: x, s,
    #          Ue_vorticity (3), Ue_max (4), Ue_fix (5),
    #          d99_vorticity (6), d99_max (7), d99_fix (8),
    #          dstar_vorticity (9), dstar_max (10), dstar_fix (11),
    #          Theta_vorticity (12), Theta_max (13), Theta_fix (14)
    Ue   = (vort=Float64.(raw[:,3]),  max=Float64.(raw[:,4]),  fix=Float64.(raw[:,5]))
    d99  = (vort=Float64.(raw[:,6]),  max=Float64.(raw[:,7]),  fix=Float64.(raw[:,8]))
    dst  = (vort=Float64.(raw[:,9]),  max=Float64.(raw[:,10]), fix=Float64.(raw[:,11]))
    Theta = (vort=Float64.(raw[:,12]), max=Float64.(raw[:,13]), fix=Float64.(raw[:,14]))

    @info "BL metrics comparison: $(length(S)) wall faces × 3 methods"

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

    colors = (vort=:royalblue, max=:firebrick, fix=:forestgreen)
    labels = (vort="vorticityIntegral", max="maxProfile", fix="fixedHeight")

    function tripanel(yvort, ymax, yfix; ylabel, title, yscale=1.0, ylog=false)
        p = plot(S, yscale .* yvort; label=labels.vort, color=colors.vort,
                 xlabel=L"S\ \mathrm{[m]}", ylabel=ylabel, title=title,
                 yscale = ylog ? :log10 : :identity, common...)
        plot!(p, S, yscale .* ymax; label=labels.max, color=colors.max, common...)
        plot!(p, S, yscale .* yfix; label=labels.fix, color=colors.fix, common...)
        return p
    end

    p1 = tripanel(Ue.vort,   Ue.max,   Ue.fix;
                  ylabel=L"U_e\ \mathrm{[m/s]}",    title="Edge velocity")
    p2 = tripanel(d99.vort,  d99.max,  d99.fix;    yscale=1e3,
                  ylabel=L"\delta_{99}\ \mathrm{[mm]}", title=L"\delta_{99}")
    p3 = tripanel(dst.vort,  dst.max,  dst.fix;    yscale=1e3,
                  ylabel=L"\delta^*\ \mathrm{[mm]}",    title=L"\delta^*")
    p4 = tripanel(Theta.vort, Theta.max, Theta.fix; yscale=1e3,
                  ylabel=L"\theta\ \mathrm{[mm]}",      title=L"\theta")

    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1400, 900))

    label = basename(case_path)
    outfile = joinpath(savedir, "BLQuantities$(label).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end
