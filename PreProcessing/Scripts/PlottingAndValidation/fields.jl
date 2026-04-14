using DelimitedFiles, LaTeXStrings

"""
    structured_reshape(x, y, f) → (xu, yu, F)

Detect structured grid from cell-center coordinates and reshape field `f`
into a matrix `F[j,i]` suitable for `heatmap`. No interpolation — pure reshaping.
"""
function structured_reshape(x::Vector{Float64}, y::Vector{Float64}, f::Vector{Float64};
                            tol_digits::Int=5)
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
    plot_fields(case_path; savedir, filename, stations, delta, strip_width,
                chord_mm, alpha_deg, x_center_mm, y_center_mm)

Parse `postProcessing/midPlane.csv` and plot velocity components (u,v,w) and pressure.
Each velocity component: scatter contour (left) + profiles at stations (right).
When airfoil geometry is provided (chord_mm != nothing), profiles are extracted
along the wall-normal direction at the given xi/c stations. Otherwise, vertical
strips at evenly spaced x-positions are used.
CSV columns: x, y, z, u, v, w, p
"""
function plot_fields(case_path::AbstractString;
                     savedir::AbstractString=case_path,
                     filename::AbstractString="midPlane.csv",
                     stations::Union{Vector{Float64}, Nothing}=nothing,
                     delta::Float64=0.03,
                     strip_width::Float64=0.005,
                     chord_mm::Union{Float64, Nothing}=nothing,
                     alpha_deg::Float64=-3.0,
                     x_center_mm::Float64=0.0,
                     y_center_mm::Float64=0.0)
    csv_path = joinpath(case_path, "postProcessing", filename)
    if !isfile(csv_path)
        @warn "midPlane.csv not found at $csv_path"
        return nothing
    end

    @info "Parsing fields from $csv_path"

    lines = readlines(csv_path)
    x_all = Float64[]; y_all = Float64[]; z_all = Float64[]
    u_all = Float64[]; v_all = Float64[]; w_all = Float64[]; p_all = Float64[]
    for line in lines[2:end]
        fields = split(line, ',')
        length(fields) == 7 || continue
        vals = tryparse.(Float64, fields)
        any(isnothing, vals) && continue
        push!(x_all, vals[1]); push!(y_all, vals[2]); push!(z_all, vals[3])
        push!(u_all, vals[4]); push!(v_all, vals[5]); push!(w_all, vals[6])
        push!(p_all, vals[7])
    end
    if isempty(x_all)
        @warn "No valid data rows in $csv_path"
        return nothing
    end

    # Filter to one z-layer (avoids NaN gaps from cyclic z-boundaries)
    z_unique = sort(unique(round.(z_all, digits=6)))
    z_target = z_unique[1]
    z_mask = abs.(z_all .- z_target) .< 1e-5
    x = x_all[z_mask]; y = y_all[z_mask]
    u = u_all[z_mask]; v = v_all[z_mask]; w = w_all[z_mask]; p = p_all[z_mask]
    @info "Parsed $(length(x_all)) cells, using $(length(x)) from z=$(round(z_target, sigdigits=4))"

    mkpath(savedir)
    label = basename(case_path)

    dx = maximum(x) - minimum(x)
    dy = maximum(y) - minimum(y)
    ar = (dx > 0 && dy > 0 && dx / dy > 5) ? :auto : :equal
    ms = length(x) > 50_000 ? 1 : 2

    use_airfoil = chord_mm !== nothing
    colors = [:royalblue, :firebrick, :forestgreen, :darkorange, :purple, :teal]

    common_sc = (
        aspect_ratio    = ar,
        colorbar        = true,
        ylabel          = L"y \ \mathrm{[m]}",
        xlabel          = L"x \ \mathrm{[m]}",
        legend          = false,
        framestyle      = :box,
        markerstrokewidth = 0,
        markersize      = ms,
        markershape     = :circle,
        tickfontsize    = 10,
        guidefontsize   = 12,
        titlefontsize   = 13,
        left_margin     = 8Plots.mm,
        bottom_margin   = 6Plots.mm,
        right_margin    = 4Plots.mm,
        top_margin      = 4Plots.mm,
        dpi             = 200,
    )

    common_prof = (
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

    # ── Profile extraction ────────────────────────────────────────────
    if use_airfoil
        prof_stations = stations !== nothing ? stations : [0.06, 0.10, 0.20, 0.30, 0.40]
        prof_ylabel   = L"\mathrm{Wall \ distance \ [m]}"
        prof_ylims    = (0.0, 0.004)

        extract_profile = (xv, yv, fv, xi_c) -> begin
            x_s, y_s, nx, ny = airfoil_surface(xi_c;
                chord_mm=chord_mm, alpha_deg=alpha_deg,
                x_center_mm=x_center_mm, y_center_mm=y_center_mm)

            raw_dn = Float64[]
            raw_dt = Float64[]
            raw_f  = Float64[]
            for i in eachindex(xv)
                rx = xv[i] - x_s
                ry = yv[i] - y_s
                d_normal  = rx * nx + ry * ny
                d_tangent = abs(rx * ny - ry * nx)
                if d_normal >= 0 && d_normal < delta && d_tangent < strip_width
                    push!(raw_dn, d_normal)
                    push!(raw_dt, d_tangent)
                    push!(raw_f,  fv[i])
                end
            end
            isempty(raw_dn) && return Float64[], Float64[]

            # Keep only the nearest tangential column
            dt_min = minimum(raw_dt)
            dt_tol = dt_min + 0.0002   # 0.2 mm tolerance
            col = raw_dt .<= dt_tol
            raw_dn = raw_dn[col]
            raw_f  = raw_f[col]

            dn_min = minimum(raw_dn)
            wall_dist = raw_dn .- dn_min
            perm = sortperm(wall_dist)
            return raw_f[perm], wall_dist[perm]
        end

        station_label = st -> latexstring("\\xi/c = $(@sprintf("%.2f", st))")

        station_x_pos = xi_c -> begin
            x_s, _, _, _ = airfoil_surface(xi_c;
                chord_mm=chord_mm, alpha_deg=alpha_deg,
                x_center_mm=x_center_mm, y_center_mm=y_center_mm)
            return x_s
        end
    else
        # Structured grid (DFP): use heatmap + column extraction
        xu, yu, F_u_s = structured_reshape(x, y, u)
        _,  _,  F_v_s = structured_reshape(x, y, v)
        _,  _,  F_w_s = structured_reshape(x, y, w)
        _,  _,  F_p_s = structured_reshape(x, y, p)
        @info "Structured grid: $(length(xu)) x $(length(yu))"

        nx_s = length(xu)
        prof_stations_idx = stations !== nothing ?
            [argmin(abs.(xu .- s)) for s in stations] :
            [1, round(Int, nx_s*0.25), round(Int, nx_s*0.50), round(Int, nx_s*0.75), nx_s]
        prof_ylabel = L"y \ \mathrm{[m]}"
        y_lo = 0.0
        y_hi = 0.25 * (maximum(yu) - minimum(yu))
        prof_ylims  = (y_lo, y_hi)

        station_label = si -> latexstring(@sprintf("x = %.3f", xu[si]), raw" \ \mathrm{m}")
        station_x_pos = si -> xu[si]
    end

    # ── Common plot options ───────────────────────────────────────────
    common_hm = (
        aspect_ratio   = ar,
        colorbar       = true,
        ylabel         = L"y \ \mathrm{[m]}",
        xlabel         = L"x \ \mathrm{[m]}",
        legend         = false,
        framestyle     = :box,
        linewidth      = 0,
        tickfontsize   = 10,
        guidefontsize  = 12,
        titlefontsize  = 13,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi            = 200,
    )

    # ── Velocity figures ──────────────────────────────────────────────
    if use_airfoil
        draw_station_lines! = (p_cont) -> begin
            for (k, st) in enumerate(prof_stations)
                x_s, y_s, nx, ny = airfoil_surface(st;
                    chord_mm=chord_mm, alpha_deg=alpha_deg,
                    x_center_mm=x_center_mm, y_center_mm=y_center_mm)
                plot!(p_cont, [x_s, x_s + delta * nx], [y_s, y_s + delta * ny];
                      color=colors[mod1(k, length(colors))],
                      linestyle=:dash, linewidth=1.2, label=false)
            end
        end

        make_velocity_figure = (fv, comp_label, comp_latex, outname) -> begin
            p_cont = scatter(x, y;
                marker_z       = fv,
                colorbar_title = latexstring(comp_label, raw" \ \mathrm{[m/s]}"),
                color          = :viridis,
                common_sc...)
            draw_station_lines!(p_cont)

            p_pr = plot(; xlabel=comp_latex, ylabel=prof_ylabel,
                ylims=prof_ylims, legend=:outerright, common_prof...)
            for (k, st) in enumerate(prof_stations)
                fprof, dprof = extract_profile(x, y, fv, st)
                isempty(fprof) && continue
                c = colors[mod1(k, length(colors))]
                plot!(p_pr, fprof, dprof;
                      label=station_label(st), color=c,
                      marker=:circle, markersize=3, markercolor=c,
                      markerstrokecolor=:black, markerstrokewidth=0.5)
            end

            fig = plot(p_cont, p_pr;
                layout = @layout([a{0.55w} b{0.45w}]), size = (1300, 450))
            outfile = joinpath(savedir, "$(outname)$(label).png")
            savefig(fig, outfile); @info "Saved: $outfile"
            return fig
        end

        fig_u = make_velocity_figure(u, "u", L"u \ \mathrm{[m/s]}", "uField")
        fig_v = make_velocity_figure(v, "v", L"v \ \mathrm{[m/s]}", "vField")
        fig_w = make_velocity_figure(w, "w", L"w \ \mathrm{[m/s]}", "wField")
        p_pres = make_velocity_figure(p, "p", L"p \ \mathrm{[m^2/s^2]}", "pressure")
    else
        make_velocity_figure = (fv_mat, comp_label, comp_latex, outname) -> begin
            p_cont = heatmap(xu, yu, fv_mat;
                colorbar_title = latexstring(comp_label, raw" \ \mathrm{[m/s]}"),
                color          = :viridis,
                common_hm...)
            for (k, si) in enumerate(prof_stations_idx)
                vline!(p_cont, [xu[si]]; color=colors[mod1(k, length(colors))],
                       linestyle=:dash, linewidth=1.2, label=false)
            end

            p_pr = plot(; xlabel=comp_latex, ylabel=prof_ylabel,
                ylims=prof_ylims, legend=:outerright, common_prof...)
            for (k, si) in enumerate(prof_stations_idx)
                profile = fv_mat[:, si]
                c = colors[mod1(k, length(colors))]
                plot!(p_pr, profile, yu;
                      label=station_label(si), color=c,
                      marker=:circle, markersize=3, markercolor=c,
                      markerstrokecolor=:black, markerstrokewidth=0.5)
            end

            fig = plot(p_cont, p_pr;
                layout = @layout([a{0.55w} b{0.45w}]), size = (1300, 450))
            outfile = joinpath(savedir, "$(outname)$(label).png")
            savefig(fig, outfile); @info "Saved: $outfile"
            return fig
        end

        fig_u = make_velocity_figure(F_u_s, "u", L"u \ \mathrm{[m/s]}", "uField")
        fig_v = make_velocity_figure(F_v_s, "v", L"v \ \mathrm{[m/s]}", "vField")
        fig_w = make_velocity_figure(F_w_s, "w", L"w \ \mathrm{[m/s]}", "wField")
        p_pres = make_velocity_figure(F_p_s, "p", L"p \ \mathrm{[m^2/s^2]}", "pressure")
    end

    return (fig_u, fig_v, fig_w, p_pres)
end
