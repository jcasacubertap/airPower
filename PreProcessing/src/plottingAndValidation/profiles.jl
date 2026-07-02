using DelimitedFiles, Printf, LaTeXStrings

# Airfoil data (ξ = x/c, η = y/c, LE → TE), upper surface only.
# Loaded from the .dat file specified in inputs.jl.
const XI_DATA, ETA_DATA = let
    dat = joinpath(ROOT, "PreProcessing", "io", "airfoilGeometryData",
                   inp.TTCP.tunnel.airfoilFile)
    load_airfoil_upper(dat)
end

# ── Global cubic spline for C2-smooth surface evaluation ──
# (must match the spline used in generateGrid.jl for mesh consistency)

struct _SplineCoeffs
    x::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
end

function _build_spline(xd::Vector{Float64}, yd::Vector{Float64})::_SplineCoeffs
    n = length(xd)
    h = [xd[i+1] - xd[i] for i in 1:n-1]
    a = copy(yd)
    lo = zeros(n); diag = zeros(n); up = zeros(n); rhs = zeros(n)
    for i in 2:n-1
        lo[i]   = h[i-1]
        diag[i] = 2.0 * (h[i-1] + h[i])
        up[i]   = h[i]
        rhs[i]  = 3.0 * ((a[i+1] - a[i]) / h[i] - (a[i] - a[i-1]) / h[i-1])
    end
    diag[1] = 1.0; up[1] = 0.0; rhs[1] = 0.0
    lo[n] = 0.0; diag[n] = 1.0; rhs[n] = 0.0
    c = zeros(n)
    for i in 2:n
        m = lo[i] / diag[i-1]; diag[i] -= m * up[i-1]; rhs[i] -= m * rhs[i-1]
    end
    c[n] = rhs[n] / diag[n]
    for i in n-1:-1:1; c[i] = (rhs[i] - up[i] * c[i+1]) / diag[i]; end
    b = zeros(n-1); d = zeros(n-1)
    for i in 1:n-1
        b[i] = (a[i+1] - a[i]) / h[i] - h[i] * (c[i+1] + 2.0 * c[i]) / 3.0
        d[i] = (c[i+1] - c[i]) / (3.0 * h[i])
    end
    return _SplineCoeffs(xd, a, b, c[1:n-1], d)
end

function _eval_spline(s::_SplineCoeffs, xq::Float64)::Float64
    n = length(s.x)
    xq ≤ s.x[1] && return s.a[1]
    xq ≥ s.x[n] && return s.a[n]
    lo, hi = 1, n - 1
    while lo < hi; mid = (lo + hi) ÷ 2; s.x[mid+1] < xq ? (lo = mid + 1) : (hi = mid); end
    dx = xq - s.x[lo]
    return s.a[lo] + dx * (s.b[lo] + dx * (s.c[lo] + dx * s.d[lo]))
end

const _ETA_SPLINE = _build_spline(XI_DATA, ETA_DATA)

"""
    airfoil_surface(xi_c; chord_mm, alpha_deg, x_center_mm, y_center_mm)

Interpolate the suction-side surface point and outward normal at chord fraction `xi_c`.
Uses a global cubic spline (C2) consistent with the mesh generator.
Returns physical coordinates [m] and unit outward normal.
"""
function airfoil_surface(xi_c::Float64;
                         chord_mm::Float64=900.0,
                         alpha_deg::Float64=-3.0,
                         x_center_mm::Float64=0.0,
                         y_center_mm::Float64=0.0)
    chord_m = chord_mm * 1e-3
    alpha   = alpha_deg * π / 180.0

    eta_local = _eval_spline(_ETA_SPLINE, xi_c)

    xb = (xi_c - 0.5) * chord_m
    yb = eta_local * chord_m
    ca, sa = cos(alpha), sin(alpha)
    x_phys =  ca * xb + sa * yb + x_center_mm * 1e-3
    y_phys = -sa * xb + ca * yb + y_center_mm * 1e-3

    # Tangent via finite difference on spline
    δ = 1e-5
    xi_lo = max(0.0, xi_c - δ)
    xi_hi = min(1.0, xi_c + δ)
    eta_lo = _eval_spline(_ETA_SPLINE, xi_lo)
    eta_hi = _eval_spline(_ETA_SPLINE, xi_hi)
    dx = (xi_hi - xi_lo) * chord_m
    dy = (eta_hi - eta_lo) * chord_m
    tx =  ca * dx + sa * dy
    ty = -sa * dx + ca * dy
    tmag = sqrt(tx^2 + ty^2)
    tx /= tmag
    ty /= tmag

    # Outward normal (suction side → rotate tangent by +90°)
    nx = -ty
    ny =  tx

    return (x_phys, y_phys, nx, ny)
end

"""
    bump_h_interp(case_path) -> Function | Nothing

Build a wall-normal bump-elevation interpolator `h(ξ/c) [m]` from the
`bumpCheck.csv` that `generateGrid.jl` writes into the case directory when the
wall bump is active (columns: xi, …, h_mm). Returns `nothing` when the bump is
disabled or the file is absent, so callers can treat the un-bumped case as
`h ≡ 0`.

The airfoil bump is parameterised in upper-surface arc length, so this CSV — the
exact samples used to deform the mesh — is the authoritative source of `h(ξ)`;
re-deriving it here would risk drifting from the mesh.

The `inp.wallModulation.enabled` guard is essential: a `bumpCheck.csv` left over
from a previous bumped run (or a cloned variant) would otherwise be picked up in
a clean case and spuriously shift the solver profile below the wall (h > 0 over
the bump support), producing a mismatch localised to the old bump's ξ-range.
"""
function bump_h_interp(case_path::AbstractString)
    inp.wallModulation.enabled || return nothing
    csv = joinpath(case_path, "bumpCheck.csv")
    isfile(csv) || return nothing
    raw = readdlm(csv, ',', Float64; skipstart=1)
    isempty(raw) && return nothing
    xi   = raw[:, 1]
    h_m  = raw[:, 7] ./ 1000.0          # h_mm → m
    xi_lo, xi_hi = xi[1], xi[end]
    return function (xi_q::Float64)
        (xi_q <= xi_lo || xi_q >= xi_hi) && return 0.0
        j = searchsortedlast(xi, xi_q)
        j = clamp(j, 1, length(xi) - 1)
        t = (xi_q - xi[j]) / (xi[j+1] - xi[j])
        return (1.0 - t) * h_m[j] + t * h_m[j+1]
    end
end

"""
    bump_support_xi(case_path) -> (xi_peak, xi_lo, xi_hi) | Nothing

Read the ξ/c support of the wall bump from `bumpCheck.csv`: the peak location
and the lower/upper chord fractions where the wall-normal elevation `h` is
non-zero. The ESN/sigmoidal shapes are truncated to exactly zero outside their
support in the CSV, so a strict `h > 0` test recovers the same extent the mesh
used. Returns `nothing` when the bump is disabled, the file is absent, or the
bump is effectively flat. (Same stale-`bumpCheck.csv` guard as `bump_h_interp`.)
"""
function bump_support_xi(case_path::AbstractString)
    inp.wallModulation.enabled || return nothing
    csv = joinpath(case_path, "bumpCheck.csv")
    isfile(csv) || return nothing
    raw = readdlm(csv, ',', Float64; skipstart=1)
    isempty(raw) && return nothing
    xi  = raw[:, 1]
    h_m = raw[:, 7] ./ 1000.0
    idx = findall(>(1e-9), h_m)
    isempty(idx) && return nothing
    return (xi[argmax(h_m)], xi[first(idx)], xi[last(idx)])
end

"""
    plot_profiles(case_path; savedir, stations, delta, filename,
                  chord_mm, alpha_deg, x_center_mm, y_center_mm)

Extract wall-normal BL profiles from midPlane.csv at specified ξ/c stations.
"""
function plot_profiles(case_path::AbstractString;
                       savedir::AbstractString=case_path,
                       stations::Vector{Float64}=[0.05, 0.10, 0.20, 0.30, 0.40],
                       delta::Float64=0.03,
                       strip_width::Float64=0.005,
                       filename::AbstractString="midPlane.csv",
                       chord_mm::Float64=900.0,
                       alpha_deg::Float64=-3.0,
                       x_center_mm::Float64=0.0,
                       y_center_mm::Float64=0.0)
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

    @info "Parsing BL profiles from $mid"
    _, raw = read_midplane(mid)
    x_all = Float64.(raw[:, 1])
    y_all = Float64.(raw[:, 2])
    u_all = Float64.(raw[:, 4])
    v_all = Float64.(raw[:, 5])
    w_all = Float64.(raw[:, 6])
    Umag_all = @. sqrt(u_all^2 + v_all^2 + w_all^2)

    # Freestream velocity estimate (99th percentile)
    sorted_U = sort(Umag_all)
    U_inf = sorted_U[round(Int, 0.99 * length(sorted_U))]

    colors = [:royalblue, :firebrick, :forestgreen, :darkorange, :purple, :teal]

    p = plot(;
        xlabel  = L"U / U_\infty",
        ylabel  = L"\mathrm{Wall \ distance \ [m]}",
        ylims   = (0.0, delta),
        size    = (680, 500),
        legend  = :outertopright,
        framestyle     = :box,
        grid           = true,
        gridalpha      = 0.3,
        gridlinewidth  = 0.5,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 9,
        titlefontsize  = 13,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi = 200,
    )

    found_any = false
    for (k, xi_c) in enumerate(stations)
        x_s, y_s, nx, ny = airfoil_surface(xi_c;
            chord_mm=chord_mm, alpha_deg=alpha_deg,
            x_center_mm=x_center_mm, y_center_mm=y_center_mm)

        raw_dn  = Float64[]
        raw_u   = Float64[]

        for i in eachindex(x_all)
            rx = x_all[i] - x_s
            ry = y_all[i] - y_s
            d_normal  = rx * nx + ry * ny
            d_tangent = abs(rx * ny - ry * nx)
            if abs(d_normal) < delta && d_tangent < strip_width
                push!(raw_dn, d_normal)
                push!(raw_u,  Umag_all[i] / U_inf)
            end
        end

        if isempty(raw_dn)
            @warn "No cells found for station xi/c = $xi_c (surface at x=$(round(x_s,digits=4)), y=$(round(y_s,digits=4)))"
            continue
        end

        # Shift so that the wall (minimum d_normal) is at y=0
        dn_min = minimum(raw_dn)
        wall_dist = raw_dn .- dn_min

        found_any = true
        perm = sortperm(wall_dist)
        c = colors[mod1(k, length(colors))]
        plot!(p, raw_u[perm], wall_dist[perm];
              label      = latexstring("\\xi/c = $(@sprintf("%.2f", xi_c))"),
              color      = c,
              linewidth  = 2,
              marker     = :circle,
              markersize = 3,
              markercolor = c,
              markerstrokewidth = 0)
    end

    if !found_any
        @warn "No profiles extracted — skipping save"
        return nothing
    end

    mkpath(savedir)
    label = basename(case_path)
    outfile = joinpath(savedir, "blProfiles$(label).png")
    savefig(p, outfile)
    @info "Saved: $outfile"
    return p
end
