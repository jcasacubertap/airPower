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
                     wm=nothing,
                     bl=nothing)
    # Accept either midPlane.csv or the binary midPlane.bin — read_midplane
    # (from blMetrics.jl) handles both; columns are x,y,z,u,v,w,p,omz.
    ppdir = joinpath(case_path, "postProcessing")
    base  = replace(filename, r"\.(csv|bin)$" => "")
    mid   = isfile(joinpath(ppdir, "$base.csv")) ? joinpath(ppdir, "$base.csv") :
            isfile(joinpath(ppdir, "$base.bin")) ? joinpath(ppdir, "$base.bin") : nothing
    if mid === nothing
        @warn "midPlane.csv/.bin not found in $ppdir"
        return nothing
    end

    @info "Parsing fields from $mid"

    _, raw = read_midplane(mid)
    if isempty(raw)
        @warn "No valid data rows in $mid"
        return nothing
    end
    x_all = Float64.(raw[:, 1]); y_all = Float64.(raw[:, 2]); z_all = Float64.(raw[:, 3])
    u_all = Float64.(raw[:, 4]); v_all = Float64.(raw[:, 5]); w_all = Float64.(raw[:, 6])
    p_all = Float64.(raw[:, 7])

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

    # ── IBL boundary-layer overlay (active when valBL=true; bl !== nothing) ──
    # Superimpose the integral-BL profile of the field named by `outname`
    # ("uField"/"vField"/"wField") at DFP streamwise position `x_st` onto the
    # given profile panel. `y_wall` shifts the BL wall-distance coordinate into
    # the panel's absolute-y frame. No-op for fields without a BL counterpart.
    bl_component = outname -> begin
        bl === nothing && return nothing
        outname == "uField" ? bl.u :
        outname == "vField" ? bl.v :
        outname == "wField" ? bl.w : nothing
    end
    # Color-matched (per-station) IBL profile — solid, thin; never labelled. A
    # single black legend entry is added once, after all station curves, by
    # bl_legend!. (DFP profiles are black lines + coloured markers, so the
    # solid coloured IBL line stays distinguishable.)
    overlay_bl! = (panel, outname, x_st, y_wall, color) -> begin
        M = bl_component(outname)
        M === nothing && return
        cidx = argmin(abs.(bl.Xgrid .- x_st))
        plot!(panel, M[:, cidx], bl.Y .+ y_wall;
              label     = false,
              color     = color,
              linewidth = 1.2)
    end
    # Single black "IBL" legend entry, placed after the station curves.
    bl_legend! = (panel, outname) -> begin
        bl_component(outname) === nothing && return
        plot!(panel, [NaN], [NaN];
              label = "IBL", color = :black, linewidth = 1.2)
    end

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

        # Wall-normal bump elevation h(ξ) [m] from bumpCheck.csv (nothing when
        # the bump is disabled). Used to anchor profiles and station markers at
        # the actual bumped wall instead of the baseline airfoil surface.
        bump_h = bump_h_interp(case_path)

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

            # Anchor at the bumped wall (subtract h(ξ)) when the bump is active;
            # otherwise fall back to the lowest-cell anchor.
            wall_dist = bump_h !== nothing ? raw_dn .- bump_h(xi_c) :
                                             raw_dn .- minimum(raw_dn)
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
                # Lift the marker base onto the bumped wall along the surface
                # normal (no-op when the bump is disabled).
                h = bump_h === nothing ? 0.0 : bump_h(st)
                x0 = x_s + h * nx; y0 = y_s + h * ny
                plot!(p_cont, [x0, x0 + delta * nx], [y0, y0 + delta * ny];
                      color=colors[mod1(k, length(colors))],
                      linestyle=:dash, linewidth=1.2, label=false)
            end
        end

        make_velocity_figure = (fv, comp_label, comp_latex, outname) -> begin
            cb_title = latexstring(comp_label, raw" \ \mathrm{[m/s]}")
            xg, yg, Fg = scatter_to_grid(x, y, fv; nx=700, ny=400, fill_iters=10)
            y_lo, y_hi = yg[1], yg[end]

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

            fig = plot(p_full, p_pr;
                layout = @layout([a{0.55w} b{0.45w}]),
                size = (1300, 450))
            outfile = joinpath(savedir, "$(outname)$(label).png")
            savefig(fig, outfile); @info "Saved: $outfile"
            return fig
        end

        fig_u = make_velocity_figure(u, "u", L"u \ \mathrm{[m/s]}", "uField")
        fig_v = make_velocity_figure(v, "v", L"v \ \mathrm{[m/s]}", "vField")
        fig_w = make_velocity_figure(w, "w", L"w \ \mathrm{[m/s]}", "wField")
        p_pres = make_velocity_figure(p, "p", L"p \ \mathrm{[m^2/s^2]}", "pressure")

        # Bump-centric profile grid (TTCP analogue of the DFP bumpProfiles panel):
        # 3 fields × 5 stations marching across the hump in ξ/c, with profiles
        # referenced to the bumped wall (extract_profile already anchors there).
        bump_sup = bump_h === nothing ? nothing : bump_support_xi(case_path)
        if bump_sup !== nothing
            xi_peak, xi_lo, xi_hi = bump_sup
            bump_w_xi   = xi_hi - xi_lo
            bump_offsets = [-1.0, -0.5, 0.0, 0.5, 1.0]
            bump_xi_stns = [clamp(xi_peak + o * bump_w_xi, 0.0, 1.0)
                            for o in bump_offsets]
            bump_titles  = [latexstring(@sprintf("\\xi/c = %.3f", st))
                            for st in bump_xi_stns]
            wd_ylabel = L"\mathrm{Wall\ distance}\ [\mathrm{m}]"

            bump_fields = [
                (u, L"u \ \mathrm{[m/s]}"),
                (v, L"v \ \mathrm{[m/s]}"),
                (w, L"w \ \mathrm{[m/s]}"),
            ]

            bump_panels = []
            for (i, (fv, comp_latex)) in enumerate(bump_fields)
                for (k, st) in enumerate(bump_xi_stns)
                    # extract_profile (use_airfoil) already returns wall distance
                    # measured from the bumped wall, so no extra h subtraction.
                    fprof, wd_prof = extract_profile(x, y, fv, st)
                    p_sub = plot(;
                        common_prof...,
                        xlabel = comp_latex,
                        ylabel = k == 1 ? wd_ylabel : "",
                        title  = i == 1 ? bump_titles[k] : "",
                        ylims  = prof_ylims,
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
                      label=station_label(st), linecolor=:black,
                      marker=:circle, markersize=3, markercolor=c,
                      markerstrokecolor=:black, markerstrokewidth=0.5)
                overlay_bl!(p_pr, outname, st, minimum(yprof), c)
            end
            bl_legend!(p_pr, outname)

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

        # ── Bump comparison figure ──────────────────────────────────────
        # Top: spanwise-velocity contour about the bump with the 5 sampling
        # stations marked. Below: u and w wall-normal profiles at those stations
        # (v omitted), one column per station, colour-keyed to the contour.
        if bump_active
            if wm.shape == :esn
                bump_xs_b, bump_xe_b = _esn_geometry_full(wm)
                xc_bump = wm.xCenter
            else
                bump_xs_b, bump_xe_b = wm.xStart, wm.xEnd
                xc_bump = (wm.xStart + wm.xEnd) / 2
            end
            bump_w_b     = bump_xe_b - bump_xs_b
            bump_offsets = [-1.0, -0.5, 0.0, 0.5, 1.0]
            bump_x_stns  = [xc_bump + o * bump_w_b for o in bump_offsets]

            # Per-station colours tie each contour line to its profile column.
            stn_cols   = [colors[mod1(k, length(colors))] for k in 1:5]
            stn_titles = [L"x_c\!-\!w_b", L"x_c\!-\!w_b/2", L"x_c",
                          L"x_c\!+\!w_b/2", L"x_c\!+\!w_b"]

            y_b_lo = minimum(y)
            wd_h   = (maximum(y) - y_b_lo) / 8          # near-wall window height
            bump_ylims = (0.0, wd_h)
            wd_ylabel  = L"\mathrm{Wall\ distance}\ [\mathrm{m}]"

            # ── Top: w contour, framed symmetrically about the bump so the 5
            #    stations sit at the 5 profile-column centres below. ──
            xg, yg, Wg = scatter_to_grid(x, y, w; nx=800, ny=320, fill_iters=20)
            xlo   = max(xg[1],   xc_bump - 1.25 * bump_w_b)
            xhi   = min(xg[end], xc_bump + 1.25 * bump_w_b)
            y_top = y_b_lo + wd_h
            p_cmap = heatmap(xg, yg, Wg;
                color          = :viridis,
                colorbar       = true,
                colorbar_title = L"w \ \mathrm{[m/s]}",
                xlims          = (xlo, xhi),
                ylims          = (yg[1], y_top),
                xlabel         = L"x \ \mathrm{[m]}",
                ylabel         = L"y \ \mathrm{[m]}",
                title          = L"\mathrm{Spanwise\ velocity\ about\ the\ wall\ bump}",
                framestyle     = :box,
                tickfontsize   = 10, guidefontsize = 12, titlefontsize = 13,
                left_margin    = 10Plots.mm, right_margin = 6Plots.mm,
                top_margin     = 5Plots.mm,  bottom_margin = 7Plots.mm,
                linewidth      = 0, dpi = 200)
            # White mask below the (bumped) wall — its upper edge draws the bump.
            xs_b = collect(range(xlo, xhi, length=600))
            ys_b = [wall_bump(xi, wm) for xi in xs_b]
            plot!(p_cmap, [xs_b; reverse(xs_b)], [ys_b; fill(yg[1], length(xs_b))];
                  seriestype=:shape, fillcolor=:white, linecolor=:white,
                  linewidth=0, label=false)
            for (k, st) in enumerate(bump_x_stns)
                (xlo <= st <= xhi) || continue
                vline!(p_cmap, [st]; color=stn_cols[k], linestyle=:dash,
                       linewidth=1.8, label=false)
                # Numbered badge linking this line to its profile column below.
                scatter!(p_cmap, [st], [y_top * 0.9]; markershape=:circle,
                         markersize=9, markercolor=stn_cols[k],
                         markerstrokecolor=:white, markerstrokewidth=1.0,
                         label=false)
                annotate!(p_cmap, st, y_top * 0.9, text(string(k), 8, :white, :center))
            end

            # ── Below: u and w profiles, one column per station ──
            # Shared x-limits per row so the stations are directly comparable.
            build_row = (fv, comp_latex, is_top) -> begin
                profs = [extract_profile(x, y, fv, st) for st in bump_x_stns]
                vmax  = maximum(pr -> isempty(pr[1]) ? 0.0 : maximum(pr[1]), profs)
                vmin  = minimum(pr -> isempty(pr[1]) ? 0.0 : minimum(pr[1]), profs)
                xr    = (min(0.0, vmin) - 0.03 * abs(vmax), 1.05 * vmax)
                row = Any[]
                for (k, st) in enumerate(bump_x_stns)
                    fprof, yprof = profs[k]
                    wd_prof = yprof .- wall_bump(st, wm)
                    p_sub = plot(; common_prof...,
                        xlabel = comp_latex,
                        ylabel = k == 1 ? wd_ylabel : "",
                        title  = is_top ? stn_titles[k] : "",
                        xlims  = xr, ylims = bump_ylims, legend = false)
                    if !isempty(fprof)
                        plot!(p_sub, fprof, wd_prof;
                              color=stn_cols[k], marker=:circle, markersize=2.6,
                              markercolor=stn_cols[k], markerstrokecolor=:black,
                              markerstrokewidth=0.3)
                    end
                    push!(row, p_sub)
                end
                return row
            end
            u_panels = build_row(u, L"u \ \mathrm{[m/s]}", true)
            w_panels = build_row(w, L"w \ \mathrm{[m/s]}", false)

            lay = @layout([ cmap{0.34h}
                            grid(2, 5) ])
            fig_bump = plot(p_cmap, u_panels..., w_panels...;
                layout = lay, size = (1650, 950), dpi = 200,
                plot_title = L"\mathrm{Base\!-\!flow\ profiles\ across\ the\ wall\ bump}",
                plot_titlefontsize = 15)
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
                      label=station_label(si), linecolor=:black,
                      marker=:circle, markersize=3, markercolor=c,
                      markerstrokecolor=:black, markerstrokewidth=0.5)
                overlay_bl!(p_pr, outname, xu[si], yu[1], c)
            end
            bl_legend!(p_pr, outname)

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

"""
    _ddy_nonuniform(y, f) → df/dy

First derivative of `f` on a (possibly non-uniform) grid `y`, both sorted by y.
Interior points use the 3-point non-uniform central stencil (2nd order); the
ends fall back to one-sided differences. Mirrors the DFP post-processor's
`diffOrder2NonUniform.m` so Julia and MATLAB agree on wall-normal derivatives.
"""
function _ddy_nonuniform(y::Vector{Float64}, f::Vector{Float64})
    n = length(y)
    d = zeros(n)
    n < 2 && return d
    d[1] = (f[2] - f[1]) / (y[2] - y[1])
    d[n] = (f[n] - f[n-1]) / (y[n] - y[n-1])
    @inbounds for i in 2:n-1
        h1 = y[i]   - y[i-1]
        h2 = y[i+1] - y[i]
        d[i] = -h2 / (h1 * (h1 + h2)) * f[i-1] +
                (h2 - h1) / (h1 * h2) * f[i]   +
                 h1 / (h2 * (h1 + h2)) * f[i+1]
    end
    return d
end

"""
    plot_dfp_first_profile(case_path; savedir, filename) → fig

DirectFlatPlate diagnostic: extract the wall-normal column at the *first cell in
x* (the inlet-most streamwise station) from `postProcessing/midPlane.{csv,bin}`
and plot the u/v/w profiles and their wall-normal derivatives ∂u/∂y, ∂v/∂y,
∂w/∂y against y. Writes `firstProfile<label>.png`. Returns the figure, or
`nothing` if the field file is missing.
"""
function plot_dfp_first_profile(case_path::AbstractString;
                                savedir::AbstractString=case_path,
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
    if isempty(raw)
        @warn "No valid data rows in $mid"
        return nothing
    end
    x_all = Float64.(raw[:, 1]); y_all = Float64.(raw[:, 2]); z_all = Float64.(raw[:, 3])
    u_all = Float64.(raw[:, 4]); v_all = Float64.(raw[:, 5]); w_all = Float64.(raw[:, 6])

    # One z-layer only (avoids duplicate columns from the cyclic span).
    z_target = sort(unique(round.(z_all, digits=6)))[1]
    zm = abs.(z_all .- z_target) .< 1e-5
    x = x_all[zm]; y = y_all[zm]; u = u_all[zm]; v = v_all[zm]; w = w_all[zm]

    # First streamwise column: snap to the smallest unique x within half the
    # spacing to the next column (robust to float noise in the cell centres).
    xu   = sort(unique(round.(x, digits=8)))
    xcol = xu[1]
    coltol = length(xu) > 1 ? 0.5 * (xu[2] - xu[1]) : 1.0
    col  = abs.(x .- xcol) .< coltol
    if !any(col)
        @warn "No cells found in the first x-column at x=$xcol"
        return nothing
    end

    yc = y[col]
    perm = sortperm(yc)
    yc = yc[perm]
    uc = u[col][perm]; vc = v[col][perm]; wc = w[col][perm]
    @info @sprintf("First-profile column at x=%.5f m: %d cells, y ∈ [%.4g, %.4g] m",
                   xcol, length(yc), yc[1], yc[end])

    dudy = _ddy_nonuniform(yc, uc)
    dvdy = _ddy_nonuniform(yc, vc)
    dwdy = _ddy_nonuniform(yc, wc)

    common = (
        framestyle     = :box,
        grid           = true,
        gridalpha      = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 9,
        titlefontsize  = 13,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        linewidth      = 2,
        dpi            = 200,
    )
    cu, cv, cw = :royalblue, :firebrick, :forestgreen

    p_vel = plot(; xlabel=L"u,\ v,\ w \ \mathrm{[m/s]}", ylabel=L"y \ \mathrm{[m]}",
                 title=latexstring(@sprintf("x = %.4f \\ \\mathrm{m}", xcol)),
                 legend=:outerright, common...)
    plot!(p_vel, uc, yc; label=L"u", color=cu, marker=:circle, markersize=2,
          markerstrokewidth=0, markercolor=cu)
    plot!(p_vel, vc, yc; label=L"v", color=cv, marker=:circle, markersize=2,
          markerstrokewidth=0, markercolor=cv)
    plot!(p_vel, wc, yc; label=L"w", color=cw, marker=:circle, markersize=2,
          markerstrokewidth=0, markercolor=cw)

    p_der = plot(; xlabel=L"\partial(u,v,w)/\partial y \ \mathrm{[1/s]}",
                 ylabel=L"y \ \mathrm{[m]}", legend=:outerright, common...)
    plot!(p_der, dudy, yc; label=L"\partial u/\partial y", color=cu,
          marker=:circle, markersize=2, markerstrokewidth=0, markercolor=cu)
    plot!(p_der, dvdy, yc; label=L"\partial v/\partial y", color=cv,
          marker=:circle, markersize=2, markerstrokewidth=0, markercolor=cv)
    plot!(p_der, dwdy, yc; label=L"\partial w/\partial y", color=cw,
          marker=:circle, markersize=2, markerstrokewidth=0, markercolor=cw)

    fig = plot(p_vel, p_der; layout=(1, 2), size=(1300, 500))

    mkpath(savedir)
    label = basename(case_path)
    outfile = joinpath(savedir, "firstProfile$(label).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end
