using DelimitedFiles, LaTeXStrings

"""
    plot_wall_quantities(case_path; savedir, chord_mm, alpha_deg,
                         x_center_mm, y_center_mm)

Plot Cp, Cf, and y+ along arc-length S from
`postProcessing/wallQuantities.csv` and `wallCoordinates.csv`.
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

    chord_m = chord_mm * 1e-3
    alpha   = alpha_deg * π / 180.0
    ca, sa  = cos(alpha), sin(alpha)
    x_ctr   = x_center_mm * 1e-3
    y_ctr   = y_center_mm * 1e-3

    raw = readdlm(csv_path, ','; skipstart=1)
    x_w  = Float64.(raw[:, 1])
    pw   = Float64.(raw[:, 2])   # p/ρ [m²/s²]
    tw   = Float64.(raw[:, 3])   # |τ_w|/ρ [m²/s²]
    yp   = Float64.(raw[:, 4])

    # Read wall coordinates for arc-length
    wc_path = joinpath(case_path, "postProcessing", "wallCoordinates.csv")
    if !isfile(wc_path)
        @warn "wallCoordinates.csv not found — using x instead of arc-length"
        perm = sortperm(x_w)
        S = x_w[perm]; pw = pw[perm]; tw = tw[perm]; yp = yp[perm]
    else
        wc = readdlm(wc_path, ','; skipstart=1)
        y_w = Float64.(wc[:, 2])

        # Sort by x/c (inverse rotation)
        xi_c = ((x_w .- x_ctr) .* ca .- (y_w .- y_ctr) .* sa) ./ chord_m .+ 0.5
        perm = sortperm(xi_c)
        x_w = x_w[perm]; y_w = y_w[perm]; pw = pw[perm]; tw = tw[perm]; yp = yp[perm]

        # Compute arc-length
        S = zeros(length(x_w))
        for i in 2:length(x_w)
            S[i] = S[i-1] + sqrt((x_w[i] - x_w[i-1])^2 + (y_w[i] - y_w[i-1])^2)
        end
    end

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

    p2 = plot(S, tw;
        xlabel = L"S \ \mathrm{[m]}",
        ylabel = L"|\tau_w|/\rho \ \mathrm{[m^2/s^2]}",
        color  = :firebrick,
        title  = "Wall shear stress",
        common...)

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
