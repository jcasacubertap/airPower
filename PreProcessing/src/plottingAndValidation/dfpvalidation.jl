using MAT, DelimitedFiles, Printf, LaTeXStrings, Statistics

"""
    plot_dfp_w_validation(case_path; savedir, gen, case_id)

w-velocity profile validation for the DirectFlatPlate case — the flat-plate
equivalent of the TTCP `plot_experimental_validation`. One subplot per
experimental station (x/c), showing the OpenFOAM w-profile (line) and the PIV
z-averaged profile (symbols), wrapped into a grid.

Reads the consolidated PIV `.mat` (a single file directly inside
`Validation/Gen{gen}/Experimental/Case{case_id}/`) holding an `output` struct —
the same schema the TTCP reader uses:
    output.xc        (1 x N)      station positions as chord percentages
    output.y         1 x N cell   each (ny x nz) wall-normal coordinate [mm]
    output.w_m_mean  1 x N cell   each (ny x nz) z-averaged spanwise velocity [m/s]

DFP-specific coordinate mapping: the flat-plate arc-length is S = xInlet + x_DFP.
Each x/c station maps to arc-length via BL.S(BL.x) (from the airfoilFlowData
`.mat`), then to x_DFP = S − xInlet, where the solver profile is sampled from a
vertical strip.
"""
function plot_dfp_w_validation(case_path::AbstractString;
                               savedir::AbstractString=case_path,
                               gen::Int=0,
                               case_id::Int=0)
    # ── Reference data for the S ↔ x mapping ──
    flow_data_dir = joinpath(ROOT, "PreProcessing", "io", "airfoilFlowData")
    fd = filter(f -> endswith(lowercase(f), ".mat"), readdir(flow_data_dir))
    if isempty(fd)
        @warn "No .mat files found in $flow_data_dir"
        return nothing
    end
    BL = matread(joinpath(flow_data_dir, fd[1]))["BL"]
    S_ref = vec(BL["S"])     # arc-length [m]
    x_ref = vec(BL["x"])     # chordwise coordinate [m]
    c_ref = BL["c"]          # chord [m]
    xInlet = inp.DFP.xInlet

    # Map x/c → x_DFP via BL.S(BL.x): S = xInlet + x_DFP.
    xc_to_xdfp = xi_c -> begin
        x_chord = xi_c * c_ref
        idx = clamp(searchsortedlast(x_ref, x_chord), 1, length(x_ref) - 1)
        t = (x_chord - x_ref[idx]) / (x_ref[idx+1] - x_ref[idx])
        S_at_xc = (1 - t) * S_ref[idx] + t * S_ref[idx+1]
        return S_at_xc - xInlet
    end

    # ── Locate the consolidated PIV .mat and read its station list ──
    val_dir = joinpath(ROOT, "PreProcessing", "io", "Validation",
                       "Gen$gen", "Experimental", "Case$case_id")
    if !isdir(val_dir)
        @warn "Validation directory not found: $val_dir"
        return nothing
    end
    mat_files = filter(f -> isfile(joinpath(val_dir, f)) && endswith(lowercase(f), ".mat"),
                       readdir(val_dir))
    if isempty(mat_files)
        @warn "No PIV .mat file in $val_dir"
        return nothing
    end
    length(mat_files) > 1 &&
        @warn "Multiple .mat files in $val_dir — using $(mat_files[1])"
    piv = matread(joinpath(val_dir, mat_files[1]))["output"]

    # xc may arrive as a Float64 matrix or a cell array of scalars / 1×1 arrays;
    # pull a plain Float64 either way.
    scalar(v) = v isa Number ? Float64(v) : Float64(first(v))
    xc_pct   = [scalar(v) for v in vec(piv["xc"])]   # chord percentages
    order    = sortperm(xc_pct)
    stations = xc_pct[order] ./ 100.0
    @info "Experimental stations (x/c): $stations"

    # ── OpenFOAM field data (accepts midPlane.csv or .bin) ──
    ppdir = joinpath(case_path, "postProcessing")
    mid   = isfile(joinpath(ppdir, "midPlane.csv")) ? joinpath(ppdir, "midPlane.csv") :
            isfile(joinpath(ppdir, "midPlane.bin")) ? joinpath(ppdir, "midPlane.bin") : nothing
    if mid === nothing
        @warn "midPlane.csv/.bin not found in $ppdir"
        return nothing
    end
    _, raw = read_midplane(mid)
    if isempty(raw)
        @warn "No valid data in $mid"
        return nothing
    end
    x_of = Float64.(raw[:, 1]); y_of = Float64.(raw[:, 2]); w_of = Float64.(raw[:, 6])

    mkpath(savedir)

    # ── Flat-plate solver profile extraction (vertical strip → wall distance) ──
    domainLength = inp.DFP.domainLength
    strip_w = 0.005 * domainLength
    extract_profile = (xv, yv, fv, x_st) -> begin
        m = abs.(xv .- x_st) .< strip_w
        any(m) || return Float64[], Float64[]
        # Snap to the nearest x-column (avoid mixing wall heights across a bump).
        dxs    = abs.(xv[m] .- x_st)
        col    = dxs .<= minimum(dxs) + 5e-5   # ~0.05 mm tolerance
        y_col  = yv[m][col]
        f_col  = fv[m][col]
        y_dist = y_col .- minimum(y_col)       # shift to wall distance
        perm   = sortperm(y_dist)
        return f_col[perm], y_dist[perm]
    end

    # ── Experimental profiles from the output cells (same as TTCP) ──
    # w_m_mean is z-averaged; reduce each (ny x nz) field to a 1-D profile by the
    # per-row valid mean. Keep the profile at its own processed y (no wall shift),
    # matching the solver which keeps its natural cell wall-distance.
    w_cells = piv["w_m_mean"]
    y_cells = piv["y"]
    exp_profiles = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    for (idx, korig) in enumerate(order)
        xi_c   = stations[idx]
        w_exp  = Float64.(w_cells[korig])          # (ny, nz) [m/s]
        y_exp  = Float64.(y_cells[korig])          # (ny, nz) [mm]
        ny_exp = size(w_exp, 1)
        w_avg  = zeros(ny_exp)
        y_dist = vec(y_exp[:, 1]) .* 1e-3          # mm → m
        for j in 1:ny_exp
            row   = w_exp[j, :]
            valid = .!isnan.(row) .& (row .!= 0)
            w_avg[j] = any(valid) ? mean(row[valid]) : NaN
        end
        keep = .!isnan.(w_avg)
        exp_profiles[xi_c] = (w_avg[keep], y_dist[keep])
        @info @sprintf("  xc=%g%%: %d points, y ∈ [%.3f, %.3f] mm",
                       xc_pct[korig], sum(keep),
                       minimum(y_dist[keep])*1e3, maximum(y_dist[keep])*1e3)
    end

    # ── Plot: one subplot per station, wrapped into a grid (same style as TTCP) ──
    N     = length(stations)
    ncols = min(N, 6)
    nrows = cld(N, ncols)
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
            ylabel = mod(k - 1, ncols) == 0 ? L"\mathrm{Wall \ distance \ [m]}" : "",
            ylims  = (0.0, 0.004),
            title  = latexstring("x/c = $(@sprintf("%.0f", st*100))", raw"\%"),
            legend = :topleft,
            common...)

        # OpenFOAM profile (line)
        x_dfp = xc_to_xdfp(st)
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
                           st, x_dfp, domainLength)
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
        layout = (nrows, ncols),
        size   = (350 * ncols, 450 * nrows),
    )

    lbl = basename(case_path)
    outfile = joinpath(savedir, "wFieldValidationExperimental$(lbl).png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
    return fig
end
