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

"""
    plot_dfp_wall_quantities(case_path; savedir, nu, filename) → fig

DirectFlatPlate analogue of `plot_wall_quantities`, in the identical 1×3 style
(wall pressure, ∂u/∂y at the wall, y⁺). The DFP case writes no
`wallQuantities.csv`, so the three quantities are reconstructed from
`postProcessing/midPlane.{csv,bin}` per streamwise x-column and plotted against
the streamwise coordinate x (the flat-plate analogue of arclength S):

  • p/ρ    — extrapolated wall value (u=v=w=0 wall point from `wallExtrapolation`)
  • ∂u/∂y  — wall gradient (u_cell / y_cell), signed; flips at separation
  • y⁺     — y1·√(|∂u/∂y|/ν), with y1 the first cell-centre wall distance
             (matches the OpenFOAM yPlus function object).
"""
function plot_dfp_wall_quantities(case_path::AbstractString;
                                  savedir::AbstractString=case_path,
                                  nu::Float64=1.456610719354608e-5,
                                  filename::AbstractString="midPlane.csv")
    ppdir = joinpath(case_path, "postProcessing")
    base  = replace(filename, r"\.(csv|bin)$" => "")
    mid   = isfile(joinpath(ppdir, "$base.csv")) ? joinpath(ppdir, "$base.csv") :
            isfile(joinpath(ppdir, "$base.bin")) ? joinpath(ppdir, "$base.bin") : nothing
    if mid === nothing
        @warn "midPlane.csv/.bin not found in $ppdir"
        return nothing
    end

    _, raw = read_midplane(mid)
    isempty(raw) && (@warn "No valid data rows in $mid"; return nothing)
    x_all = Float64.(raw[:, 1]); y_all = Float64.(raw[:, 2]); z_all = Float64.(raw[:, 3])
    u_all = Float64.(raw[:, 4]); p_all = Float64.(raw[:, 7])

    # One z-layer only (avoids duplicate columns from the cyclic span).
    z_target = sort(unique(round.(z_all, digits=6)))[1]
    zm = abs.(z_all .- z_target) .< 1e-5
    x = x_all[zm]; y = y_all[zm]; u = u_all[zm]; p = p_all[zm]

    # Group cells into streamwise columns by (rounded) x.
    xr = round.(x, digits=8)
    cols = Dict{Float64, Vector{Int}}()
    for i in eachindex(xr)
        push!(get!(cols, xr[i], Int[]), i)
    end

    xw = Float64[]; pw = Float64[]; dudy = Float64[]; yp = Float64[]
    for (xc, idx) in cols
        length(idx) < 2 && continue
        o  = sortperm(y[idx])
        yy = y[idx][o]; uu = u[idx][o]; pp = p[idx][o]
        # Wall gradient from the wall point (yy[1]≈0, uu[1]≈0) to the first cell.
        g  = (uu[2] - uu[1]) / (yy[2] - yy[1])
        y1 = yy[2] - yy[1]                       # first cell-centre wall distance
        utau = sqrt(nu * abs(g))
        push!(xw, xc); push!(pw, pp[1]); push!(dudy, g); push!(yp, utau * y1 / nu)
    end
    isempty(xw) && (@warn "No wall columns reconstructed from $mid"; return nothing)

    # Order by streamwise coordinate for monotone plots.
    perm = sortperm(xw)
    xw = xw[perm]; pw = pw[perm]; dudy = dudy[perm]; yp = yp[perm]
    @info "DFP wall quantities: $(length(xw)) columns"

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

    p1 = plot(xw, pw;
        xlabel = L"x \ \mathrm{[m]}",
        ylabel = L"p/\rho \ \mathrm{[m^2/s^2]}",
        color  = :royalblue,
        title  = "Wall pressure",
        common...)

    p2 = plot(xw, dudy;
        xlabel = L"x \ \mathrm{[m]}",
        ylabel = L"\partial u/\partial y \ \mathrm{[1/s]}",
        color  = :firebrick,
        title  = L"\partial u/\partial y \ \mathrm{at\ wall\ (sign\ flips\ at\ separation)}",
        common...)
    Plots.hline!(p2, [0]; color=:black, linestyle=:dash, linewidth=1, label="")

    p3 = plot(xw, yp;
        xlabel = L"x \ \mathrm{[m]}",
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
