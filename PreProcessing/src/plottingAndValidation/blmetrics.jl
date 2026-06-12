using DelimitedFiles, LaTeXStrings, MAT

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
    plot_ue_xc(case_path; savedir, ref_mat)

Plot the edge velocity U_e versus x/c — the CFD (from `BLQuantities.csv`,
configured method) and, when `ref_mat` is given, the reference U_e from the
flow-data `.mat` (`BL["x"]/BL["c"]`, `BL["Ue"]`). x/c is common to both, so the
comparison is direct (no arclength-origin reconciliation). Saves to
`<savedir>/UeComparison<caseLabel>.png`.
"""
function plot_ue_xc(case_path::AbstractString;
                    savedir::AbstractString=case_path,
                    ref_mat::Union{AbstractString, Nothing}=nothing)
    csv_path = joinpath(case_path, "postProcessing", "BLQuantities.csv")
    if !isfile(csv_path)
        @warn "BLQuantities.csv not found at $csv_path"
        return nothing
    end
    raw = readdlm(csv_path, ','; skipstart=1)   # x, s, xi, Ue, d99, dstar, Theta
    xi  = Float64.(raw[:, 3])
    Ue  = Float64.(raw[:, 4])
    mkpath(savedir)

    p = plot(xi, Ue;
        label="CFD", color=:royalblue, linewidth=2,
        xlabel=L"x/c", ylabel=L"U_e\ \mathrm{[m/s]}", title="Edge velocity",
        framestyle=:box, grid=true, gridalpha=0.3,
        tickfontsize=10, guidefontsize=12, titlefontsize=13,
        left_margin=8Plots.mm, bottom_margin=6Plots.mm,
        legend=:bottomright, dpi=200, size=(800, 500))

    if ref_mat !== nothing && isfile(ref_mat)
        BL   = matread(ref_mat)["BL"]
        cref = BL["c"]; cref = cref isa Number ? Float64(cref) : Float64(first(cref))
        xr   = vec(Float64.(BL["x"])) ./ cref
        ur   = vec(Float64.(BL["Ue"]))
        o    = sortperm(xr); xr = xr[o]; ur = ur[o]
        keep = (xr .>= minimum(xi)) .& (xr .<= maximum(xi))     # clip to CFD range
        if any(keep)
            plot!(p, xr[keep], ur[keep];
                  label="Reference (.mat)", color=:black, linestyle=:dash,
                  linewidth=2, marker=:circle, markersize=3, markerstrokewidth=0)
        end
    end

    outfile = joinpath(savedir, "UeComparison$(basename(case_path)).png")
    savefig(p, outfile)
    @info "Saved: $outfile"
    return p
end
