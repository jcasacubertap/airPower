using MAT, DelimitedFiles, Printf, LaTeXStrings

"""
    plot_experimental_validation(case_path; savedir, gen, chord_mm, alpha_deg,
                                 x_center_mm, y_center_mm, delta, strip_width)

One subplot per experimental station (x/c), each showing OpenFOAM w-profile
(line) and PIV w_m z-averaged profile (symbols).  Stations are discovered from
`Validation/Gen{gen}/Experimental/Case0/xc{NN}/stats_raw.mat`.
"""
function plot_experimental_validation(case_path::AbstractString;
                                      savedir::AbstractString=case_path,
                                      gen::Int=0,
                                      chord_mm::Float64=900.0,
                                      alpha_deg::Float64=-3.0,
                                      x_center_mm::Float64=0.0,
                                      y_center_mm::Float64=0.0,
                                      delta::Float64=0.010,
                                      strip_width::Float64=0.005)
    # ── Discover experimental stations ──
    val_dir = joinpath(ROOT, "PreProcessing", "InputOutput", "Validation",
                       "Gen$gen", "Experimental", "Case0")
    if !isdir(val_dir)
        @warn "Validation directory not found: $val_dir"
        return nothing
    end

    station_dirs = filter(d -> startswith(d, "xc") && isdir(joinpath(val_dir, d)),
                          readdir(val_dir))
    if isempty(station_dirs)
        @warn "No xc* station folders in $val_dir"
        return nothing
    end

    # Parse x/c fractions from folder names (e.g. "xc15" → 0.15)
    stations = Float64[]
    for sd in station_dirs
        num = tryparse(Int, replace(sd, "xc" => ""))
        num === nothing && continue
        push!(stations, num / 100.0)
    end
    sort!(stations)
    @info "Experimental stations (x/c): $stations"

    # ── Read OpenFOAM field data ──
    csv_path = joinpath(case_path, "postProcessing", "midPlane.csv")
    if !isfile(csv_path)
        @warn "midPlane.csv not found at $csv_path"
        return nothing
    end

    lines = readlines(csv_path)
    x_of = Float64[]; y_of = Float64[]; w_of = Float64[]
    for line in lines[2:end]
        fields = split(line, ',')
        length(fields) == 7 || continue
        vals = tryparse.(Float64, fields)
        any(isnothing, vals) && continue
        push!(x_of, vals[1]); push!(y_of, vals[2]); push!(w_of, vals[6])
    end
    if isempty(x_of)
        @warn "No valid data in $csv_path"
        return nothing
    end

    mkpath(savedir)

    # ── Profile extraction (same logic as fields.jl) ──
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

        dt_min = minimum(raw_dt)
        dt_tol = dt_min + 0.0002
        col = raw_dt .<= dt_tol
        raw_dn = raw_dn[col]
        raw_f  = raw_f[col]

        dn_min = minimum(raw_dn)
        wall_dist = raw_dn .- dn_min
        perm = sortperm(wall_dist)
        return raw_f[perm], wall_dist[perm]
    end

    # ── Read experimental profiles ──
    exp_profiles = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    for xi_c in stations
        folder = joinpath(val_dir, "xc$(Int(round(xi_c * 100)))")
        mat_file = joinpath(folder, "stats_raw.mat")
        if !isfile(mat_file)
            @warn "Missing $mat_file"
            continue
        end

        piv = matread(mat_file)["stats_piv_raw"]["c1"]
        w_exp = Float64.(piv["w_m"])   # (ny, nz), spanwise velocity [m/s]
        y_exp = Float64.(piv["y"])     # (ny, nz), mm

        # Average w_m along z (dim 2), ignoring zeros/NaNs
        ny_exp = size(w_exp, 1)
        w_avg  = zeros(ny_exp)
        y_avg  = y_exp[:, 1] .* 1e-3   # mm → m

        for j in 1:ny_exp
            row = w_exp[j, :]
            valid = .!isnan.(row) .& (row .!= 0)
            if any(valid)
                w_avg[j] = mean(row[valid])
            else
                w_avg[j] = NaN
            end
        end

        # Locate the wall: y where w is minimum (no-slip → w ≈ 0)
        valid = .!isnan.(w_avg)
        wall_idx = argmin(w_avg[valid])
        y_wall = y_avg[valid][wall_idx]

        # Shift so wall = 0, keep only above wall
        y_shifted = y_avg .- y_wall
        mask = valid .& (y_shifted .>= 0)
        exp_profiles[xi_c] = (w_avg[mask], y_shifted[mask])
        @info @sprintf("  xc%d: %d points, wall at y = %.3f mm", Int(round(xi_c*100)), sum(mask), y_wall*1e3)
    end

    # ── Plot: one subplot per station ──
    N = length(stations)
    common = (
        framestyle     = :box,
        grid           = true,
        gridalpha      = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 8,
        titlefontsize  = 11,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi            = 200,
    )

    subplots = []
    for (k, st) in enumerate(stations)
        p = plot(;
            xlabel = L"w \ \mathrm{[m/s]}",
            ylabel = k == 1 ? L"\mathrm{Wall \ distance \ [m]}" : "",
            ylims  = (0.0, 0.004),
            title  = latexstring("x/c = $(@sprintf("%.0f", st*100))", raw"\%"),
            legend = :topleft,
            common...)

        # OpenFOAM profile (line)
        wprof, dprof = extract_profile(x_of, y_of, w_of, st)
        if !isempty(wprof)
            plot!(p, wprof, dprof;
                label="Numerical", color=:black, linewidth=2)
        end

        # Experimental profile (symbols)
        if haskey(exp_profiles, st)
            w_e, y_e = exp_profiles[st]
            scatter!(p, w_e, y_e;
                label="Experimental", color=:firebrick,
                marker=:circle, markersize=2, markerstrokewidth=0)
        end

        push!(subplots, p)
    end

    fig = plot(subplots...;
        layout = (1, N),
        size   = (350 * N, 450),
    )

    lbl = basename(case_path)
    outfile = joinpath(savedir, "wFieldValidationExperimental$(lbl).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end
