"""
    parse_patch_face_centers(case_path, patch_name) → (x, y)

Extract face centres of the named patch from the OpenFOAM mesh.
Returns mid-plane face centres sorted by x.
"""
function parse_patch_face_centers(case_path::AbstractString, patch_name::AbstractString)
    poly = joinpath(case_path, "constant", "polyMesh")

    # ── Read points ──
    pts_lines = readlines(joinpath(poly, "points"))
    npts = 0
    data_start = 0
    for (i, line) in enumerate(pts_lines)
        m = match(r"^(\d+)$", strip(line))
        if m !== nothing && npts == 0
            npts = parse(Int, m.captures[1])
        end
        if strip(line) == "("
            data_start = i + 1
            break
        end
    end

    px = Vector{Float64}(undef, npts)
    py = Vector{Float64}(undef, npts)
    pz = Vector{Float64}(undef, npts)
    for j in 1:npts
        m = match(r"\(\s*([^ ]+)\s+([^ ]+)\s+([^ )]+)\s*\)", pts_lines[data_start + j - 1])
        px[j] = parse(Float64, m.captures[1])
        py[j] = parse(Float64, m.captures[2])
        pz[j] = parse(Float64, m.captures[3])
    end

    # ── Find patch in boundary ──
    bnd_text = read(joinpath(poly, "boundary"), String)
    bm = match(Regex("$(patch_name)\\s*\\{[^}]*nFaces\\s+(\\d+)[^}]*startFace\\s+(\\d+)"), bnd_text)
    if bm === nothing
        error("Patch '$patch_name' not found in boundary file")
    end
    nFaces    = parse(Int, bm.captures[1])
    startFace = parse(Int, bm.captures[2])

    # ── Read faces ──
    faces_lines = readlines(joinpath(poly, "faces"))
    f_data_start = 0
    for (i, line) in enumerate(faces_lines)
        if strip(line) == "("
            f_data_start = i + 1
            break
        end
    end

    fx = Float64[]
    fy = Float64[]
    fz = Float64[]
    for k in 1:nFaces
        line = faces_lines[f_data_start + startFace + k - 1]
        m = match(r"\d+\(([^)]+)\)", line)
        verts = parse.(Int, split(m.captures[1]))
        cx = sum(px[v+1] for v in verts) / length(verts)
        cy = sum(py[v+1] for v in verts) / length(verts)
        cz = sum(pz[v+1] for v in verts) / length(verts)
        push!(fx, cx); push!(fy, cy); push!(fz, cz)
    end

    # Filter to one z-layer (mid-plane)
    z_vals = sort(unique(round.(fz, digits=8)))
    if length(z_vals) > 1
        mask = abs.(fz .- z_vals[1]) .< 0.5 * abs(z_vals[2] - z_vals[1])
        fx = fx[mask]; fy = fy[mask]
    end

    perm = sortperm(fx)
    return fx[perm], fy[perm]
end

"""
    plot_wall_geometry(case_path; savedir, chord_mm, alpha_deg,
                       x_center_mm, y_center_mm)

Plot three airfoil coordinate sets for validation:
  1. Original .dat file (transformed to physical coordinates)
  2. Mesh airfoil patch face centres
  3. Exported wallCoordinates.csv from postProcessing
"""
function plot_wall_geometry(case_path::AbstractString;
                            savedir::AbstractString=case_path,
                            chord_mm::Float64=900.0,
                            alpha_deg::Float64=-3.0,
                            x_center_mm::Float64=0.0,
                            y_center_mm::Float64=0.0)
    mkpath(savedir)

    # ── 1. Original airfoil from .dat file (full contour, markers only) ──
    chord_m = chord_mm * 1e-3
    alpha   = alpha_deg * π / 180.0
    ca, sa  = cos(alpha), sin(alpha)
    x_ctr   = x_center_mm * 1e-3
    y_ctr   = y_center_mm * 1e-3

    dat_path = joinpath(ROOT, "PreProcessing", "InputOutput", "AirfoilGeometryData",
                        inp.TTCP.tunnel.airfoilFile)
    x_orig = Float64[]
    y_orig = Float64[]
    for line in eachline(dat_path)
        line = strip(line)
        isempty(line) && continue
        parts = split(line)
        length(parts) >= 2 || continue
        xi = tryparse(Float64, parts[1])
        eta = tryparse(Float64, parts[2])
        (xi === nothing || eta === nothing) && continue
        xb = (xi - 0.5) * chord_m
        yb = eta * chord_m
        push!(x_orig, x_ctr + ca * xb + sa * yb)
        push!(y_orig, y_ctr - sa * xb + ca * yb)
    end
    @info "Original airfoil: $(length(x_orig)) points"

    # ── 2. Mesh airfoil patch face centres ──
    x_mesh, y_mesh = Float64[], Float64[]
    try
        x_mesh, y_mesh = parse_patch_face_centers(case_path, "airfoil")
        @info "Mesh airfoil: $(length(x_mesh)) face centres"
    catch e
        @warn "Could not parse mesh airfoil patch: $e"
    end

    # ── 3. Exported wall coordinates ──
    x_out, y_out = Float64[], Float64[]
    wc_path = joinpath(case_path, "postProcessing", "wallCoordinates.csv")
    if isfile(wc_path)
        raw = readdlm(wc_path, ','; skipstart=1)
        x_out = Float64.(raw[:, 1])
        y_out = Float64.(raw[:, 2])
        perm = sortperm(x_out)
        x_out = x_out[perm]
        y_out = y_out[perm]
        @info "Output wall: $(length(x_out)) points"
    else
        @warn "wallCoordinates.csv not found at $wc_path"
    end

    # ── LE zoom limits (10% chord region around LE) ──
    chord_m = chord_mm * 1e-3
    x_le, y_le, _, _ = airfoil_surface(0.0;
        chord_mm=chord_mm, alpha_deg=alpha_deg,
        x_center_mm=x_center_mm, y_center_mm=y_center_mm)
    le_span = 0.05 * chord_m
    xlims_zoom = (x_le - 0.3 * le_span, x_le + le_span)
    ylims_zoom = (y_le - 0.6 * le_span, y_le + 0.6 * le_span)

    # ── Common plot options ──
    common = (
        aspect_ratio = :equal,
        xlabel     = L"x \ \mathrm{[m]}",
        ylabel     = L"y \ \mathrm{[m]}",
        framestyle = :box,
        grid       = true,
        gridalpha  = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 9,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi        = 200,
    )

    # Helper to add mesh + output data to a subplot
    function add_data!(pl)
        if !isempty(x_mesh)
            scatter!(pl, x_mesh, y_mesh;
                label      = "Mesh (patch faces)",
                color      = :royalblue,
                markersize = 2,
                markerstrokewidth = 0,
            )
        end
        if !isempty(x_out)
            scatter!(pl, x_out, y_out;
                label      = "Output (wallCoordinates)",
                color      = :firebrick,
                markersize = 3,
                markershape = :xcross,
                markerstrokewidth = 1,
            )
        end
    end

    # ── Full view ──
    p1 = plot(x_orig, y_orig;
        label = "Original (.dat)", color = :black,
        seriestype = :scatter, markersize = 3, markerstrokewidth = 0,
        legend = :topright, title = "Full airfoil",
        titlefontsize = 11, common...)
    add_data!(p1)

    # ── LE zoom (same x-span as full view, no aspect_ratio constraint) ──
    p2 = plot(x_orig, y_orig;
        label = "Original (.dat)", color = :black,
        seriestype = :scatter, markersize = 3, markerstrokewidth = 0,
        legend = false, title = "LE zoom (5% chord)",
        titlefontsize = 11,
        xlims = xlims_zoom, ylims = ylims_zoom,
        aspect_ratio = :auto,
        xlabel     = L"x \ \mathrm{[m]}",
        ylabel     = L"y \ \mathrm{[m]}",
        framestyle = :box,
        grid       = true,
        gridalpha  = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 9,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi        = 200,
    )
    add_data!(p2)

    fig = plot(p1, p2;
        layout = grid(2, 1),
        size   = (900, 800),
    )

    label = basename(case_path)
    outfile = joinpath(savedir, "wallGeometry$(label).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end

"""
    plot_airfoil_bump_check(savedir; csv_path)

Read the (xi, x_orig, y_orig, x_bumped, y_bumped, s, h) samples written by
generateGrid.jl when the wall bump is active and produce a three-panel figure:
    • top:    full upper-surface contour with both curves overlaid
    • middle: zoom around the bump region
    • bottom: wall-normal elevation h(s) of the bumped surface relative to
              the original, with peak height and width annotated

Returns the figure, or `nothing` if `bumpCheck.csv` is missing (e.g. bump
was disabled when meshAirfoil ran).
"""
function plot_airfoil_bump_check(; savedir::AbstractString,
                                   csv_path::AbstractString)
    if !isfile(csv_path)
        @warn "bumpCheck.csv not found — skipping airfoil bump check plot" path=csv_path
        return nothing
    end
    mkpath(savedir)

    raw = readdlm(csv_path, ','; skipstart=1)
    xi_arr = Float64.(raw[:, 1])             # chord fraction
    x_orig = Float64.(raw[:, 2]) ./ 1000.0   # mm → m
    y_orig = Float64.(raw[:, 3]) ./ 1000.0
    x_bump = Float64.(raw[:, 4]) ./ 1000.0
    y_bump = Float64.(raw[:, 5]) ./ 1000.0
    s_arr  = Float64.(raw[:, 6])             # [mm]
    h_arr  = Float64.(raw[:, 7])             # [mm]

    # Locate the bump support by point-wise displacement
    disp = hypot.(x_bump .- x_orig, y_bump .- y_orig)
    idx_bump = findall(d -> d > 1e-9, disp)

    common = (
        aspect_ratio = :equal,
        xlabel       = L"x \ \mathrm{[m]}",
        ylabel       = L"y \ \mathrm{[m]}",
        framestyle   = :box,
        grid         = true,
        gridalpha    = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 9,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi            = 200,
    )

    p_full = plot(x_orig, y_orig;
        label = "Original upper surface", color = :black, linewidth = 3,
        legend = :topright, title = "Full upper surface", titlefontsize = 11,
        common...)
    plot!(p_full, x_bump, y_bump;
        label = "Bumped upper surface", color = :firebrick, linewidth = 1)

    if isempty(idx_bump)
        @warn "bumpCheck.csv contains no perturbed samples — bump effectively zero"
        fig = p_full
    else
        # Zoom: extend 5× the bump span on each side, but no smaller than 10×
        # peak displacement in y so the bump shape is clearly visible.
        i_lo = max(1, first(idx_bump))
        i_hi = min(length(disp), last(idx_bump))
        x_b_lo, x_b_hi = extrema(x_orig[i_lo:i_hi])
        bump_x_span = x_b_hi - x_b_lo
        peak_disp   = maximum(disp[i_lo:i_hi])
        pad_x = 0.5 * bump_x_span
        pad_y = max(0.5 * bump_x_span, 10 * peak_disp)
        y_b_mid = 0.5 * (minimum(y_orig[i_lo:i_hi]) + maximum(y_orig[i_lo:i_hi]))
        xlims_zoom = (x_b_lo - pad_x, x_b_hi + pad_x)
        ylims_zoom = (y_b_mid - pad_y, y_b_mid + pad_y)

        p_zoom = plot(x_orig, y_orig;
            label = "Original upper surface", color = :black, linewidth = 3,
            legend = :topright, title = "Bump region zoom", titlefontsize = 11,
            xlims = xlims_zoom, ylims = ylims_zoom,
            aspect_ratio = :auto,
            xlabel = L"x \ \mathrm{[m]}", ylabel = L"y \ \mathrm{[m]}",
            framestyle = :box, grid = true, gridalpha = 0.3,
            tickfontsize = 10, guidefontsize = 12, legendfontsize = 9,
            left_margin = 8Plots.mm, bottom_margin = 6Plots.mm,
            right_margin = 4Plots.mm, top_margin = 4Plots.mm, dpi = 200)
        plot!(p_zoom, x_bump, y_bump;
            label = "Bumped upper surface", color = :firebrick, linewidth = 1)

        # Third panel: wall-normal elevation h(s) of the bumped surface
        # relative to the original.  Convert s [mm] → s [m] for display.
        i_peak_arr = argmax(abs.(h_arr))
        peak_h_mm  = h_arr[i_peak_arr]              # signed peak (handles depressions)
        peak_xi    = xi_arr[i_peak_arr]             # chord fraction at the peak
        bump_width_mm = s_arr[last(idx_bump)] - s_arr[first(idx_bump)]
        s_pad = 0.5 * bump_width_mm
        s_lo  = max(0.0, s_arr[first(idx_bump)] - s_pad)
        s_hi  = s_arr[last(idx_bump)] + s_pad
        annot_text = @sprintf("h_max = %.3f mm   width = %.2f mm   (x/c) at peak = %.4f",
                              peak_h_mm, bump_width_mm, peak_xi)

        p_h = plot(s_arr ./ 1000.0, h_arr;
            label = false, color = :firebrick, linewidth = 1.5,
            title = "Wall-normal elevation:  $annot_text",
            titlefontsize = 11,
            xlims = (s_lo / 1000.0, s_hi / 1000.0),
            xlabel = L"s \ \mathrm{[m]}",
            ylabel = L"h \ \mathrm{[mm]}",
            framestyle = :box, grid = true, gridalpha = 0.3,
            aspect_ratio = :auto,
            tickfontsize = 10, guidefontsize = 12,
            left_margin = 8Plots.mm, bottom_margin = 6Plots.mm,
            right_margin = 4Plots.mm, top_margin = 4Plots.mm, dpi = 200)
        hline!(p_h, [0.0]; color = :black, linewidth = 3, label = false)

        # yTol-truncated bounds (dashed black) — first/last s where |h| > yTol
        yTol_mm = inp.wallModulation.yTol * 1000.0
        idx_tol = findall(x -> abs(x) > yTol_mm, h_arr)
        if !isempty(idx_tol)
            s_tol_lo_m = s_arr[first(idx_tol)] / 1000.0
            s_tol_hi_m = s_arr[last(idx_tol)]  / 1000.0
            vline!(p_h, [s_tol_lo_m, s_tol_hi_m];
                color = :black, linestyle = :dash, linewidth = 1, label = false)
        end

        # Real arc-length center (dotted black) — peak of |h|
        i_peak     = argmax(abs.(h_arr))
        s_center_m = s_arr[i_peak] / 1000.0
        vline!(p_h, [s_center_m];
            color = :black, linestyle = :dot, linewidth = 1, label = false)

        fig = plot(p_full, p_zoom, p_h; layout = grid(3, 1), size = (900, 1100))
    end

    outfile = joinpath(savedir, "airfoilWallModulationCheck.png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end
