using DelimitedFiles, LaTeXStrings

"""
    structured_reshape(x, y, f) → (xu, yu, F)

Detect structured grid from cell-center coordinates and reshape field `f`
into a matrix `F[j,i]` suitable for `heatmap`. No interpolation — pure reshaping.
"""
function structured_reshape(x::Vector{Float64}, y::Vector{Float64}, f::Vector{Float64};
                            tol_digits::Int=8)
    xr = round.(x, digits=tol_digits)
    yr = round.(y, digits=tol_digits)
    xu = sort(unique(xr))
    yu = sort(unique(yr))

    xi_map = Dict(v => i for (i, v) in enumerate(xu))
    yi_map = Dict(v => j for (j, v) in enumerate(yu))

    F = fill(NaN, length(yu), length(xu))
    for k in eachindex(x)
        i = xi_map[xr[k]]
        j = yi_map[yr[k]]
        F[j, i] = f[k]
    end
    return xu, yu, F
end

"""
    plot_fields(case_path; savedir, filename)

Parse `postProcessing/midPlane.csv` and plot velocity components (u,v,w) and pressure.
Each velocity component: contour (left) + wall-normal profiles at 5 stations (right).
CSV columns: x, y, z, u, v, w, p
"""
function plot_fields(case_path::AbstractString;
                     savedir::AbstractString=case_path,
                     filename::AbstractString="midPlane.csv")
    csv_path = joinpath(case_path, "postProcessing", filename)
    if !isfile(csv_path)
        @warn "midPlane.csv not found at $csv_path"
        return nothing
    end

    @info "Parsing fields from $csv_path"

    lines = readlines(csv_path)
    x = Float64[]; y = Float64[]; u = Float64[]; v = Float64[]
    w = Float64[]; p = Float64[]
    for line in lines[2:end]
        fields = split(line, ',')
        length(fields) == 7 || continue
        vals = tryparse.(Float64, fields)
        any(isnothing, vals) && continue
        push!(x, vals[1]); push!(y, vals[2])
        push!(u, vals[4]); push!(v, vals[5]); push!(w, vals[6])
        push!(p, vals[7])
    end
    if isempty(x)
        @warn "No valid data rows in $csv_path"
        return nothing
    end
    @info "Parsed $(length(x)) cells"

    mkpath(savedir)
    label = basename(case_path)

    dx = maximum(x) - minimum(x)
    dy = maximum(y) - minimum(y)
    ar = (dx > 0 && dy > 0 && dx / dy > 5) ? :auto : :equal

    xu, yu, F_u = structured_reshape(x, y, u)
    _,  _,  F_v = structured_reshape(x, y, v)
    _,  _,  F_w = structured_reshape(x, y, w)
    _,  _,  F_p = structured_reshape(x, y, p)
    @info "Structured grid: $(length(xu)) x $(length(yu))"

    # 5 profile stations: inlet, outlet, and 3 equi-spaced in between
    nx = length(xu)
    station_idx = [1,
                   round(Int, nx * 0.25),
                   round(Int, nx * 0.50),
                   round(Int, nx * 0.75),
                   nx]

    colors = [:royalblue, :firebrick, :forestgreen, :darkorange, :purple]

    common_hm = (
        aspect_ratio   = ar,
        colorbar       = true,
        ylabel         = L"y \ \mathrm{[m]}",
        xlabel         = L"x \ \mathrm{[m]}",
        legend         = false,
        framestyle     = :box,
        tickfontsize   = 10,
        guidefontsize  = 12,
        titlefontsize  = 13,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi            = 200,
    )

    # y-range for profile plots: from y=0 up to quarter of the domain height
    y_lo = 0.0
    y_hi = 0.25 * (maximum(yu) - minimum(yu))

    common_prof = (
        ylabel         = L"y \ \mathrm{[m]}",
        ylims          = (y_lo, y_hi),
        framestyle     = :box,
        grid           = true,
        gridalpha      = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 7,
        titlefontsize  = 13,
        left_margin    = 5Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        linewidth      = 2,
        dpi            = 200,
    )

    # Helper: generate contour + profiles figure for one velocity component
    function make_velocity_figure(F, comp_label, comp_latex, outname)
        p_hm = heatmap(xu, yu, F;
            colorbar_title = latexstring(comp_label, raw" \ \mathrm{[m/s]}"),
            color          = :viridis,
            common_hm...)

        # Add vertical dashed lines at profile stations
        for (k, si) in enumerate(station_idx)
            vline!(p_hm, [xu[si]]; color=colors[k], linestyle=:dash, linewidth=1.2, label=false)
        end

        p_pr = plot(;
            xlabel = comp_latex,
            legend = :outerright,
            common_prof...)

        for (k, si) in enumerate(station_idx)
            profile = F[:, si]
            lbl = latexstring(@sprintf("x = %.3f", xu[si]), raw" \ \mathrm{m}")
            plot!(p_pr, profile, yu;
                  label=lbl, color=colors[k])
        end

        fig = plot(p_hm, p_pr;
            layout = @layout([a{0.55w} b{0.45w}]),
            size   = (1300, 450),
        )

        outfile = joinpath(savedir, "$(outname)$(label).png")
        savefig(fig, outfile)
        @info "Saved: $outfile"
        return fig
    end

    fig_u = make_velocity_figure(F_u, "u", L"u \ \mathrm{[m/s]}", "uField")
    fig_v = make_velocity_figure(F_v, "v", L"v \ \mathrm{[m/s]}", "vField")
    fig_w = make_velocity_figure(F_w, "w", L"w \ \mathrm{[m/s]}", "wField")

    # --- Pressure (standalone) ---
    p_pres = heatmap(xu, yu, F_p;
        colorbar_title = L"p \ \mathrm{[Pa]}",
        color          = :RdBu,
        clims          = (-maximum(abs.(p)), maximum(abs.(p))),
        size           = (950, 420),
        common_hm...)
    outfile_p = joinpath(savedir, "pressure$(label).png")
    savefig(p_pres, outfile_p)
    @info "Saved: $outfile_p"

    return (fig_u, fig_v, fig_w, p_pres)
end
