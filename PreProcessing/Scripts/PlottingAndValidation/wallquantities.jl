using DelimitedFiles, LaTeXStrings

"""
    plot_wall_quantities(case_path; savedir, chord_mm, alpha_deg,
                         x_center_mm, y_center_mm)

Plot Cp, Cf, and y+ along arc-length S from `postProcessing/wallQuantities.csv`.
The arclength `s` is now a column of that file (written by writeMidPlane), so
this function no longer needs `wallCoordinates.csv`.

The chord_mm/alpha/center args are kept for backward compatibility but unused.
"""
function plot_wall_quantities(case_path::AbstractString;
                              savedir::AbstractString=case_path,
                              chord_mm::Float64=900.0,
                              alpha_deg::Float64=-3.0,
                              x_center_mm::Float64=0.0,
                              y_center_mm::Float64=0.0)
    csv_path = joinpath(case_path, "postProcessing", "wallQuantities.csv")
    if !isfile(csv_path)
        @warn "wallQuantities.csv not found at $csv_path"
        return nothing
    end

    raw = readdlm(csv_path, ','; skipstart=1)
    x_w  = Float64.(raw[:, 1])
    S    = Float64.(raw[:, 2])   # arclength along upper surface from xi=0 [m]
    pw   = Float64.(raw[:, 3])   # p/ρ [m²/s²]
    dudy = Float64.(raw[:, 4])   # ∂u/∂y at the wall [1/s] (signed; flips at separation)
    yp   = Float64.(raw[:, 5])

    # Order by arclength for monotone plots
    perm = sortperm(S)
    S = S[perm]; pw = pw[perm]; dudy = dudy[perm]; yp = yp[perm]

    @info "Wall quantities: $(length(S)) points"

    mkpath(savedir)

    common = (
        framestyle     = :box,
        grid           = true,
        gridalpha      = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        titlefontsize  = 13,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        legend         = false,
        linewidth      = 2,
        dpi            = 200,
    )

    p1 = plot(S, pw;
        xlabel = L"S \ \mathrm{[m]}",
        ylabel = L"p/\rho \ \mathrm{[m^2/s^2]}",
        color  = :royalblue,
        title  = "Wall pressure",
        common...)

    p2 = plot(S, dudy;
        xlabel = L"S \ \mathrm{[m]}",
        ylabel = L"\partial u/\partial y \ \mathrm{[1/s]}",
        color  = :firebrick,
        title  = L"\partial u/\partial y \ \mathrm{at\ wall\ (sign\ flips\ at\ separation)}",
        common...)
    Plots.hline!(p2, [0]; color=:black, linestyle=:dash, linewidth=1, label="")

    p3 = plot(S, yp;
        xlabel = L"S \ \mathrm{[m]}",
        ylabel = L"y^+",
        color  = :forestgreen,
        title  = L"y^+ \ \mathrm{distribution}",
        common...)

    fig = plot(p1, p2, p3;
        layout = (1, 3),
        size   = (1500, 400),
    )

    label = basename(case_path)
    outfile = joinpath(savedir, "wallQuantities$(label).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end
