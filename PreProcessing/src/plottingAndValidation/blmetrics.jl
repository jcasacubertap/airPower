using DelimitedFiles, LaTeXStrings

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
    plot_bl_metrics_comparison(case_path; savedir)

Read `postProcessing/BLQuantities_compare.csv` and produce a 2×2 figure
overlaying the five U_e methods (vorticityIntegralTrapezoidal,
vorticityIntegralMidpoint, maxProfile, fixedHeight, pressureBernoulli)
on each panel:

    ┌─────────────────┬─────────────────┐
    │     U_e(s)      │    δ99(s)       │
    ├─────────────────┼─────────────────┤
    │     δ*(s)       │     θ(s)        │
    └─────────────────┴─────────────────┘

Saves to `<savedir>/BLQuantitiesCompare<caseLabel>.png` and returns the
figure.
"""
function plot_bl_metrics_comparison(case_path::AbstractString;
                                    savedir::AbstractString=case_path)
    csv_path = joinpath(case_path, "postProcessing", "BLQuantities_compare.csv")
    if !isfile(csv_path)
        @warn "BLQuantities_compare.csv not found at $csv_path"
        return nothing
    end

    raw = readdlm(csv_path, ','; skipstart=1)
    S   = Float64.(raw[:, 2])
    # Columns: x, s,
    #          Ue_vorticityTrap (3), Ue_vorticityMid (4), Ue_max (5), Ue_fix (6), Ue_bern (7),
    #          d99_*   (8..12), dstar_* (13..17), Theta_* (18..22)
    Ue    = (vortT=Float64.(raw[:,3]), vortM=Float64.(raw[:,4]),
             max=Float64.(raw[:,5]),   fix=Float64.(raw[:,6]),  bern=Float64.(raw[:,7]))
    d99   = (vortT=Float64.(raw[:,8]), vortM=Float64.(raw[:,9]),
             max=Float64.(raw[:,10]),  fix=Float64.(raw[:,11]), bern=Float64.(raw[:,12]))
    dst   = (vortT=Float64.(raw[:,13]),vortM=Float64.(raw[:,14]),
             max=Float64.(raw[:,15]),  fix=Float64.(raw[:,16]), bern=Float64.(raw[:,17]))
    Theta = (vortT=Float64.(raw[:,18]),vortM=Float64.(raw[:,19]),
             max=Float64.(raw[:,20]),  fix=Float64.(raw[:,21]), bern=Float64.(raw[:,22]))

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

    colors = (vortT=:royalblue, vortM=:darkorange, max=:firebrick,
              fix=:forestgreen, bern=:purple)
    labels = (vortT="vorticityIntegralTrapezoidal", vortM="vorticityIntegralMidpoint",
              max="maxProfile", fix="fixedHeight", bern="pressureBernoulli")

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
