using DelimitedFiles, LaTeXStrings

"""
    scatter_to_grid(x, y, f; nx, ny, fill_iters) → (xg, yg, F)

Bin scattered cell-center data to a regular nx×ny grid, averaging samples that
fall in the same bin. Empty bins are then filled by iterative neighbor-averaging
(`fill_iters` passes), which closes free-stream gaps where cells are larger than
a bin. `fill_iters` is capped so genuine no-data regions (e.g. airfoil interior)
stay NaN and render as the plot background.
"""
function scatter_to_grid(x::Vector{Float64}, y::Vector{Float64}, f::Vector{Float64};
                         nx::Int=600, ny::Int=300, fill_iters::Int=15)
    x_lo, x_hi = extrema(x)
    y_lo, y_hi = extrema(y)
    dx_b = (x_hi - x_lo) / max(nx - 1, 1)
    dy_b = (y_hi - y_lo) / max(ny - 1, 1)
    xg = collect(range(x_lo, x_hi, length=nx))
    yg = collect(range(y_lo, y_hi, length=ny))

    F_sum = zeros(ny, nx)
    F_cnt = zeros(Int, ny, nx)
    @inbounds for k in eachindex(x)
        i = clamp(round(Int, (x[k] - x_lo) / dx_b) + 1, 1, nx)
        j = clamp(round(Int, (y[k] - y_lo) / dy_b) + 1, 1, ny)
        F_sum[j, i] += f[k]
        F_cnt[j, i] += 1
    end

    F = fill(NaN, ny, nx)
    @inbounds for k in eachindex(F)
        F_cnt[k] > 0 && (F[k] = F_sum[k] / F_cnt[k])
    end

    for _ in 1:fill_iters
        any(isnan, F) || break
        F_new = copy(F)
        changed = false
        @inbounds for j in 1:ny, i in 1:nx
            isnan(F[j, i]) || continue
            s = 0.0; c = 0
            for dj in -1:1, di in -1:1
                (di == 0 && dj == 0) && continue
                ii = i + di; jj = j + dj
                (1 <= ii <= nx && 1 <= jj <= ny) || continue
                v = F[jj, ii]
                isnan(v) && continue
                s += v; c += 1
            end
            if c > 0
                F_new[j, i] = s / c
                changed = true
            end
        end
        F = F_new
        changed || break
    end
    return xg, yg, F
end

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
                     y_center_mm::Float64=0.0,
                     wm=nothing)
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

    use_airfoil = chord_mm !== nothing

    # Detect if mesh is structured (for heatmap) or deformed (use scatter)
    use_scatter = use_airfoil  # airfoil always uses scatter
    if !use_scatter
        # Check if structured reshape would have NaNs (deformed mesh, e.g., bump)
        xu_test, yu_test, F_test = structured_reshape(x, y, u)
        nan_frac = sum(isnan.(F_test)) / length(F_test)
        if nan_frac > 0.01
            use_scatter = true
            @info "Deformed mesh detected ($(round(nan_frac*100, digits=1))% NaN) — using scatter plots"
        end
    end
    colors = [:royalblue, :firebrick, :forestgreen, :darkorange, :purple, :teal]

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
    elseif use_scatter
        # Deformed DFP grid (bump enabled): scatter + vertical strip profiles
        x_lo_v, x_hi_v = extrema(x)
        prof_stations = stations !== nothing ? stations : collect(range(x_lo_v, x_hi_v, length=5))
        strip_w_vert = 0.005 * dx
        prof_ylabel = L"y \ \mathrm{[m]}"
        prof_ylims  = (minimum(y), minimum(y) + 0.25 * dy)

        extract_profile = (xv, yv, fv, x_st) -> begin
            in_strip = abs.(xv .- x_st) .< strip_w_vert
            !any(in_strip) && return Float64[], Float64[]
            # Snap to the single x-column nearest x_st (avoid interleaving
            # multiple columns when strip_w_vert spans more than one).
            dxs = abs.(xv .- x_st)
            dx_min = minimum(dxs[in_strip])
            col_tol = max(1e-9, 0.01 * strip_w_vert)
            mask = in_strip .& (dxs .<= dx_min + col_tol)
            perm = sortperm(yv[mask])
            return fv[mask][perm], yv[mask][perm]
        end

        station_label = st -> latexstring(@sprintf("x = %.3f", st), raw" \ \mathrm{m}")
        station_x_pos = st -> st
    else
        # Structured grid (flat DFP): heatmap + column extraction
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
            cb_title = latexstring(comp_label, raw" \ \mathrm{[m/s]}")
            xg, yg, Fg = scatter_to_grid(x, y, fv; nx=700, ny=400, fill_iters=10)
            y_lo, y_hi = yg[1], yg[end]
            y_zoom_hi = y_lo + (y_hi - y_lo) / 8

            make_panel = (ylims_p; force_auto_ar=false, show_cb_title=true) -> begin
                panel = heatmap(xg, yg, Fg;
                    common_hm...,
                    colorbar_title = show_cb_title ? cb_title : "",
                    color          = :viridis,
                    xlims          = (xg[1], xg[end]),
                    ylims          = ylims_p,
                    aspect_ratio   = force_auto_ar ? :auto : ar)
                draw_station_lines!(panel)
                return panel
            end

            p_full = make_panel((y_lo, y_hi))
            p_zoom = make_panel((y_lo, y_zoom_hi); force_auto_ar=true, show_cb_title=false)

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

            fig = plot(p_full, p_zoom, p_pr;
                layout = @layout([grid(2,1){0.55w} b{0.45w}]),
                size = (1300, 450))
            outfile = joinpath(savedir, "$(outname)$(label).png")
            savefig(fig, outfile); @info "Saved: $outfile"
            return fig
        end

        fig_u = make_velocity_figure(u, "u", L"u \ \mathrm{[m/s]}", "uField")
        fig_v = make_velocity_figure(v, "v", L"v \ \mathrm{[m/s]}", "vField")
        fig_w = make_velocity_figure(w, "w", L"w \ \mathrm{[m/s]}", "wField")
        p_pres = make_velocity_figure(p, "p", L"p \ \mathrm{[m^2/s^2]}", "pressure")
    elseif use_scatter
        # Deformed DFP (bump): gridded heatmap + vertical strip profiles
        bump_active = wm !== nothing && hasproperty(wm, :enabled) && wm.enabled
        make_velocity_figure = (fv, comp_label, comp_latex, outname) -> begin
            cb_title = latexstring(comp_label, raw" \ \mathrm{[m/s]}")
            xg, yg, Fg = scatter_to_grid(x, y, fv; nx=700, ny=300, fill_iters=20)
            y_lo, y_hi = yg[1], yg[end]
            y_zoom_hi = y_lo + (y_hi - y_lo) / 8

            make_panel = (ylims_p; show_cb_title=true) -> begin
                panel = heatmap(xg, yg, Fg;
                    common_hm...,
                    colorbar_title = show_cb_title ? cb_title : "",
                    color          = :viridis,
                    xlims          = (xg[1], xg[end]),
                    ylims          = ylims_p)
                if bump_active
                    xs_b = collect(range(xg[1], xg[end], length=600))
                    ys_b = [wall_bump(xi, wm) for xi in xs_b]
                    plot!(panel,
                          [xs_b; reverse(xs_b)],
                          [ys_b; fill(yg[1], length(xs_b))];
                          seriestype=:shape, fillcolor=:white,
                          linecolor=:white, linewidth=0, label=false)
                end
                for (k, st) in enumerate(prof_stations)
                    vline!(panel, [st]; color=colors[mod1(k, length(colors))],
                           linestyle=:dash, linewidth=1.2, label=false)
                end
                return panel
            end

            p_full = make_panel((y_lo, y_hi))
            p_zoom = make_panel((y_lo, y_zoom_hi); show_cb_title=false)

            p_pr = plot(; xlabel=comp_latex, ylabel=prof_ylabel,
                ylims=prof_ylims, legend=:outerright, common_prof...)
            for (k, st) in enumerate(prof_stations)
                fprof, yprof = extract_profile(x, y, fv, st)
                isempty(fprof) && continue
                c = colors[mod1(k, length(colors))]
                plot!(p_pr, fprof, yprof;
                      label=station_label(st), color=c,
                      marker=:circle, markersize=3, markercolor=c,
                      markerstrokecolor=:black, markerstrokewidth=0.5)
            end

            fig = plot(p_full, p_zoom, p_pr;
                layout = @layout([grid(2,1){0.55w} b{0.45w}]),
                size = (1300, 450))
            outfile = joinpath(savedir, "$(outname)$(label).png")
            savefig(fig, outfile); @info "Saved: $outfile"
            return fig
        end

        fig_u = make_velocity_figure(u, "u", L"u \ \mathrm{[m/s]}", "uField")
        fig_v = make_velocity_figure(v, "v", L"v \ \mathrm{[m/s]}", "vField")
        fig_w = make_velocity_figure(w, "w", L"w \ \mathrm{[m/s]}", "wField")
        p_pres = make_velocity_figure(p, "p", L"p \ \mathrm{[m^2/s^2]}", "pressure")

        # Bump-centric profile grid (4 fields × 5 stations spaced by bump width)
        if bump_active
            if wm.shape == :esn
                bump_xs_b, bump_xe_b = _esn_geometry_full(wm)
                xc_bump = wm.xCenter
            else
                bump_xs_b, bump_xe_b = wm.xStart, wm.xEnd
                xc_bump = (wm.xStart + wm.xEnd) / 2
            end
            bump_w_b = bump_xe_b - bump_xs_b
            bump_offsets = [-1.0, -0.5, 0.0, 0.5, 1.0]
            bump_x_stns  = [xc_bump + o * bump_w_b for o in bump_offsets]
            bump_titles  = [latexstring(@sprintf("x = %.4f \\ \\mathrm{m}", st))
                            for st in bump_x_stns]

            y_b_lo = minimum(y); y_b_hi = maximum(y)
            bump_ylims = (0.0, (y_b_hi - y_b_lo) / 8)
            wd_ylabel  = L"\mathrm{Wall\ distance}\ [\mathrm{m}]"

            bump_fields = [
                (u, L"u \ \mathrm{[m/s]}"),
                (v, L"v \ \mathrm{[m/s]}"),
                (w, L"w \ \mathrm{[m/s]}"),
            ]

            bump_panels = []
            for (i, (fv, comp_latex)) in enumerate(bump_fields)
                for (k, st) in enumerate(bump_x_stns)
                    fprof, yprof = extract_profile(x, y, fv, st)
                    y_wall_st = wall_bump(st, wm)
                    wd_prof = yprof .- y_wall_st
                    p_sub = plot(;
                        common_prof...,
                        xlabel = comp_latex,
                        ylabel = k == 1 ? wd_ylabel : "",
                        title  = i == 1 ? bump_titles[k] : "",
                        ylims  = bump_ylims,
                        legend = false)
                    if !isempty(fprof)
                        plot!(p_sub, fprof, wd_prof;
                              color             = :royalblue,
                              marker            = :circle,
                              markersize        = 3,
                              markercolor       = :royalblue,
                              markerstrokecolor = :black,
                              markerstrokewidth = 0.5)
                    end
                    push!(bump_panels, p_sub)
                end
            end

            fig_bump = plot(bump_panels...;
                layout = (3, 5), size = (1500, 800), dpi = 200)
            outfile = joinpath(savedir, "bumpProfiles$(label).png")
            savefig(fig_bump, outfile)
            @info "Saved: $outfile"
        end
    else
        make_velocity_figure = (fv_mat, comp_label, comp_latex, outname) -> begin
            cb_title = latexstring(comp_label, raw" \ \mathrm{[m/s]}")
            y_lo, y_hi = minimum(yu), maximum(yu)
            y_zoom_hi = y_lo + (y_hi - y_lo) / 8

            make_panel = (ylims_p; show_cb_title=true) -> begin
                panel = heatmap(xu, yu, fv_mat;
                    common_hm...,
                    colorbar_title = show_cb_title ? cb_title : "",
                    color          = :viridis,
                    ylims          = ylims_p)
                for (k, si) in enumerate(prof_stations_idx)
                    vline!(panel, [xu[si]]; color=colors[mod1(k, length(colors))],
                           linestyle=:dash, linewidth=1.2, label=false)
                end
                return panel
            end

            p_full = make_panel((y_lo, y_hi))
            p_zoom = make_panel((y_lo, y_zoom_hi); show_cb_title=false)

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

            fig = plot(p_full, p_zoom, p_pr;
                layout = @layout([grid(2,1){0.55w} b{0.45w}]),
                size = (1300, 450))
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
