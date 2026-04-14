using MAT, DelimitedFiles, Printf, LaTeXStrings, Statistics

"""
    plot_dfp_w_validation(case_path; savedir, gen)

w-velocity profile validation for the DirectFlatPlate case.
One subplot per experimental station (x/c), showing OpenFOAM w-profile (line)
and PIV w_m z-averaged profile (symbols).

Coordinate mapping: the flat plate arc-length is S = xInlet + x_OpenFOAM.
The x/c station maps to arc-length via BL.S(BL.x), then to x_DFP = S - xInlet.
"""
function plot_dfp_w_validation(case_path::AbstractString;
                               savedir::AbstractString=case_path,
                               gen::Int=0)
    # ── Load reference data for S ↔ x mapping ──
    flow_data_dir = joinpath(ROOT, "PreProcessing", "InputOutput", "AirfoilFlowData")
    mat_files = filter(f -> endswith(f, ".mat"), readdir(flow_data_dir))
    if isempty(mat_files)
        @warn "No .mat files found in $flow_data_dir"
        return nothing
    end
    BL = matread(joinpath(flow_data_dir, mat_files[1]))["BL"]
    S_ref = vec(BL["S"])     # arc-length [m]
    x_ref = vec(BL["x"])     # chordwise coordinate [m]
    c_ref = BL["c"]          # chord [m]
    xInlet = inp.DFP.xInlet

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

    stations = Float64[]
    for sd in station_dirs
        num = tryparse(Int, replace(sd, "xc" => ""))
        num === nothing && continue
        push!(stations, num / 100.0)
    end
    sort!(stations)
    @info "Experimental stations (x/c): $stations"

    # ── Map x/c → x_DFP via BL.S(BL.x) ──
    # Domain starts at x=0, physical arc-length S = xInlet + x_OpenFOAM
    function xc_to_xdfp(xi_c)
        x_chord = xi_c * c_ref
        # Interpolate S at this x_chord
        idx = searchsortedlast(x_ref, x_chord)
        idx = clamp(idx, 1, length(x_ref) - 1)
        t = (x_chord - x_ref[idx]) / (x_ref[idx+1] - x_ref[idx])
        S_at_xc = (1 - t) * S_ref[idx] + t * S_ref[idx+1]
        return S_at_xc - xInlet
    end

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

    # ── Flat plate profile extraction (vertical strip) ──
    domainLength = inp.DFP.domainLength
    strip_w = 0.005 * domainLength

    extract_profile = (xv, yv, fv, x_st) -> begin
        mask = abs.(xv .- x_st) .< strip_w
        !any(mask) && return Float64[], Float64[]
        perm = sortperm(yv[mask])
        return fv[mask][perm], yv[mask][perm]
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
        w_exp = Float64.(piv["w_m"])
        y_exp = Float64.(piv["y"])

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

        valid = .!isnan.(w_avg)
        wall_idx = argmin(w_avg[valid])
        y_wall = y_avg[valid][wall_idx]

        y_shifted = y_avg .- y_wall
        mask = valid .& (y_shifted .>= 0)
        exp_profiles[xi_c] = (w_avg[mask], y_shifted[mask])
        @info @sprintf("  xc%d: %d points, wall at y = %.3f mm",
                        Int(round(xi_c*100)), sum(mask), y_wall*1e3)
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

    domainHeight = inp.DFP.domainHeight
    subplots = []
    for (k, xi_c) in enumerate(stations)
        x_dfp = xc_to_xdfp(xi_c)
        @info @sprintf("  x/c=%.2f → x_DFP=%.4f m", xi_c, x_dfp)

        p = plot(;
            xlabel = L"w \ \mathrm{[m/s]}",
            ylabel = k == 1 ? L"y \ \mathrm{[m]}" : "",
            ylims  = (0.0, 0.004),
            title  = latexstring("x/c = $(@sprintf("%.0f", xi_c*100))", raw"\%"),
            legend = :topleft,
            common...)

        # OpenFOAM profile
        if 0 < x_dfp < domainLength
            wprof, yprof = extract_profile(x_of, y_of, w_of, x_dfp)
            if !isempty(wprof)
                plot!(p, wprof, yprof;
                    label="Solver", color=:black, linewidth=2,
                    marker=:circle, markersize=3, markercolor=:black,
                    markerstrokecolor=:black, markerstrokewidth=0.5)
            end
        else
            @warn @sprintf("  x/c=%.2f maps to x_DFP=%.4f — outside domain [0, %.3f]",
                           xi_c, x_dfp, domainLength)
        end

        # Experimental profile
        if haskey(exp_profiles, xi_c)
            w_e, y_e = exp_profiles[xi_c]
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
