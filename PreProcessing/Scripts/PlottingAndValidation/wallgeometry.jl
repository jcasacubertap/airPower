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
