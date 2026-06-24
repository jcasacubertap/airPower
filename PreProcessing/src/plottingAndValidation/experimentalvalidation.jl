using MAT, DelimitedFiles, Printf, LaTeXStrings

"""
    plot_experimental_validation(case_path; savedir, gen, chord_mm, alpha_deg,
                                 x_center_mm, y_center_mm, delta, strip_width)

One subplot per experimental station (x/c), each showing the OpenFOAM w-profile
(line) and the PIV w_m_mean profile (symbols).  Stations and profiles are read
from the single `.mat` inside
`Validation/Gen{gen}/Experimental/Case{case_id}/`, holding an `output` struct:
    output.xc        (1 x N)      station positions as chord percentages
    output.y         1 x N cell   each (ny x nz) wall-normal coordinate [mm]
    output.w_m_mean  1 x N cell   each (ny x nz) z-averaged spanwise velocity [m/s]
"""
function plot_experimental_validation(case_path::AbstractString;
                                      savedir::AbstractString=case_path,
                                      gen::Int=0,
                                      case_id::Int=0,
                                      chord_mm::Float64=900.0,
                                      alpha_deg::Float64=-3.0,
                                      x_center_mm::Float64=0.0,
                                      y_center_mm::Float64=0.0,
                                      delta::Float64=0.010,
                                      strip_width::Float64=0.005)
    # ── Locate the PIV .mat and read its station list ──
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

    # xc may arrive as a Float64 matrix (Case0) or as a cell array of scalars /
    # 1×1 arrays (Case1, Matrix{Any}); pull a plain Float64 either way.
    scalar(v) = v isa Number ? Float64(v) : Float64(first(v))
    xc_pct   = [scalar(v) for v in vec(piv["xc"])]   # chord percentages, e.g. 20 → x/c=0.20
    order    = sortperm(xc_pct)
    stations = xc_pct[order] ./ 100.0
    @info "Experimental stations (x/c): $stations"

    # ── Read OpenFOAM field data ──
    # Accept either midPlane.csv or the binary midPlane.bin — read_midplane
    # (from blMetrics.jl) handles both; columns are x,y,z,u,v,w,p,omz.
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

    # Wall-normal bump elevation h(ξ) [m] (nothing when the bump is disabled).
    # The PIV wall-normal coordinate is anchored at the *real* (bumped) wall, so
    # the solver's distance — measured from the baseline airfoil surface — is
    # shifted down by h(ξ) to share that reference.
    bump_h = bump_h_interp(case_path)

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

        # Reference the wall-normal distance to the bumped wall (h above the
        # baseline surface), matching the PIV anchor. No-op when h ≡ 0.
        bump_h !== nothing && (raw_dn = raw_dn .- bump_h(xi_c))

        # natural wall-normal distance of each cell from the (bumped) wall
        perm = sortperm(raw_dn)
        return raw_f[perm], raw_dn[perm]
    end

    # ── Read experimental profiles from the .mat cells ──
    # w_m_mean is already z-averaged (constant across z); reduce each (ny x nz)
    # field to a 1-D profile by the per-row valid mean (a no-op for the
    # pre-averaged data, robust to z dropouts). The profile keeps its own
    # processed y (output.y) — no wall anchoring, matching the CFD which keeps
    # its natural cell wall-distance.
    w_cells = piv["w_m_mean"]
    y_cells = piv["y"]
    exp_profiles = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    for (idx, korig) in enumerate(order)
        xi_c  = stations[idx]
        w_exp = Float64.(w_cells[korig])   # (ny, nz), spanwise velocity [m/s]
        y_exp = Float64.(y_cells[korig])   # (ny, nz), mm (constant across z)

        ny_exp = size(w_exp, 1)
        w_avg  = zeros(ny_exp)
        y_dist = vec(y_exp[:, 1]) .* 1e-3   # mm → m, the PIV-processed wall-normal coordinate

        for j in 1:ny_exp
            row   = w_exp[j, :]
            valid = .!isnan.(row) .& (row .!= 0)
            w_avg[j] = any(valid) ? mean(row[valid]) : NaN
        end

        # keep the profile at its own processed y — no wall search, no shift
        keep = .!isnan.(w_avg)
        exp_profiles[xi_c] = (w_avg[keep], y_dist[keep])
        @info @sprintf("  xc=%g%%: %d points, y ∈ [%.3f, %.3f] mm",
                       xc_pct[korig], sum(keep), minimum(y_dist[keep])*1e3, maximum(y_dist[keep])*1e3)
    end

    # ── Plot: one subplot per station (wrapped into a grid if many) ──
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
        wprof, dprof = extract_profile(x_of, y_of, w_of, st)
        if !isempty(wprof)
            plot!(p, wprof, dprof;
                label="Solver", color=:black, linewidth=2,
                marker=:circle, markersize=3, markercolor=:black,
                markerstrokecolor=:black, markerstrokewidth=0.5)
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
