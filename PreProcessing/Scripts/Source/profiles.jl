using DelimitedFiles, Printf

# M3J airfoil data (ξ = x/c, η = y/c, LE → TE), 65 points on suction side.
# Copied from generateGrid.jl.
const XI_DATA = Float64[
    0.0,      0.001313, 0.003646, 0.007137, 0.011779,
    0.017562, 0.024472, 0.032492, 0.041604, 0.051787,
    0.063016, 0.075266, 0.088508, 0.10271,  0.11784,
    0.13387,  0.15074,  0.16844,  0.18691,  0.20611,
    0.22599,  0.24652,  0.26764,  0.2893,   0.31145,
    0.33405,  0.35703,  0.38034,  0.40394,  0.42776,
    0.45174,  0.47584,  0.5,      0.52416,  0.54826,
    0.57224,  0.59606,  0.61966,  0.64297,  0.66595,
    0.68855,  0.7107,   0.73236,  0.75348,  0.77401,
    0.79389,  0.81309,  0.83156,  0.84926,  0.86613,
    0.88216,  0.89729,  0.91149,  0.92473,  0.93698,
    0.94821,  0.9584,   0.96751,  0.97553,  0.98244,
    0.98822,  0.99286,  0.99635,  0.99869,  1.0
]

const ETA_DATA = Float64[
    0.0,      0.005524, 0.009305, 0.013007, 0.016561,
    0.019935, 0.023174, 0.026436, 0.029709, 0.033033,
    0.036384, 0.039752, 0.043132, 0.046513, 0.049875,
    0.053194, 0.056452, 0.059635, 0.062735, 0.065716,
    0.068591, 0.071357, 0.073995, 0.076469, 0.078799,
    0.080952, 0.082919, 0.084694, 0.086235, 0.087554,
    0.088594, 0.089378, 0.089851, 0.09002,  0.089814,
    0.08923,  0.088209, 0.086598, 0.084245, 0.080977,
    0.077081, 0.07263,  0.067808, 0.062635, 0.057286,
    0.051831, 0.046364, 0.040993, 0.035807, 0.03088,
    0.026254, 0.021944, 0.017969, 0.014393, 0.011285,
    0.00868,  0.006551, 0.004917, 0.003622, 0.002599,
    0.001744, 0.001056, 0.00054,  0.000194, 0.0
]

"""
    interp_airfoil(xi_c) → (x_phys, y_phys, nx, ny)

Interpolate the suction-side surface point and outward normal at chord fraction `xi_c`.
Returns physical coordinates [m] and unit outward normal.
Reads chord, alpha, center from AirfoilLECase/system/inputDomain.
"""
function airfoil_surface(xi_c::Float64;
                         chord_mm::Float64=900.0,
                         alpha_deg::Float64=-3.0,
                         x_center_mm::Float64=0.0,
                         y_center_mm::Float64=0.0)
    chord_m = chord_mm * 1e-3
    alpha   = alpha_deg * π / 180.0

    # Linear interpolation in M3J data
    n = length(XI_DATA)
    # Find bracketing interval
    idx = searchsortedlast(XI_DATA, xi_c)
    idx = clamp(idx, 1, n - 1)
    t = (xi_c - XI_DATA[idx]) / (XI_DATA[idx+1] - XI_DATA[idx])

    xi_local  = (1 - t) * XI_DATA[idx]  + t * XI_DATA[idx+1]
    eta_local = (1 - t) * ETA_DATA[idx] + t * ETA_DATA[idx+1]

    # Body-frame → physical: x_body = (ξ - 0.5)*chord, y_body = η*chord, then rotate by alpha
    xb = (xi_local - 0.5) * chord_m
    yb = eta_local * chord_m
    ca, sa = cos(alpha), sin(alpha)
    x_phys = ca * xb - sa * yb + x_center_mm * 1e-3
    y_phys = sa * xb + ca * yb + y_center_mm * 1e-3

    # Tangent via finite difference on M3J data
    if idx < n - 1
        xi_next  = XI_DATA[idx+1]
        eta_next = ETA_DATA[idx+1]
    else
        xi_next  = XI_DATA[idx]
        eta_next = ETA_DATA[idx]
    end
    if idx > 1
        xi_prev  = XI_DATA[idx]
        eta_prev = ETA_DATA[idx]
    else
        xi_prev  = XI_DATA[idx+1]
        eta_prev = ETA_DATA[idx+1]
    end
    dx = (xi_next - xi_prev) * chord_m
    dy = (eta_next - eta_prev) * chord_m
    # Rotate tangent
    tx = ca * dx - sa * dy
    ty = sa * dx + ca * dy
    tmag = sqrt(tx^2 + ty^2)
    tx /= tmag
    ty /= tmag

    # Outward normal (suction side → rotate tangent by +90°)
    nx = -ty
    ny =  tx

    return (x_phys, y_phys, nx, ny)
end

"""
    plot_profiles(case_path; savedir, stations, delta, filename,
                  chord_mm, alpha_deg, x_center_mm, y_center_mm)

Extract wall-normal BL profiles from midPlane.csv at specified ξ/c stations.
`stations` are fractional chord positions (0 = LE, 1 = TE).
`delta` is the max wall-normal distance to include [m].
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
    csv_path = joinpath(case_path, "postProcessing", filename)
    if !isfile(csv_path)
        @warn "midPlane.csv not found at $csv_path"
        return nothing
    end

    @info "Parsing BL profiles from $csv_path"
    raw = readdlm(csv_path, ','; skipstart=1)
    x_all = Float64.(raw[:, 1])
    y_all = Float64.(raw[:, 2])
    u_all = Float64.(raw[:, 4])
    v_all = Float64.(raw[:, 5])
    w_all = Float64.(raw[:, 6])
    Umag_all = @. sqrt(u_all^2 + v_all^2 + w_all^2)

    # Freestream velocity estimate (99th percentile to avoid outliers)
    sorted_U = sort(Umag_all)
    U_inf = sorted_U[round(Int, 0.99 * length(sorted_U))]

    p = plot(; xlabel="U / U∞", ylabel="Wall distance [m]",
             title="BL profiles — $(basename(case_path))",
             size=(800, 600), legend=:topright, grid=true)

    found_any = false
    for xi_c in stations
        x_s, y_s, nx, ny = airfoil_surface(xi_c;
            chord_mm=chord_mm, alpha_deg=alpha_deg,
            x_center_mm=x_center_mm, y_center_mm=y_center_mm)

        # Collect all cells in the strip (allow large normal range to handle LE interpolation offset)
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
            @warn "No cells found for station ξ/c = $xi_c (surface at x=$(round(x_s,digits=4)), y=$(round(y_s,digits=4)))"
            continue
        end

        # Shift so that the wall (minimum d_normal) is at y=0
        dn_min = minimum(raw_dn)
        wall_dist = raw_dn .- dn_min

        found_any = true
        perm = sortperm(wall_dist)
        plot!(p, raw_u[perm], wall_dist[perm];
              label=@sprintf("ξ/c = %.2f", xi_c), linewidth=1.5, marker=:circle, markersize=2)
    end

    if !found_any
        @warn "No profiles extracted — skipping save"
        return nothing
    end

    mkpath(savedir)
    label = basename(case_path)
    outfile = joinpath(savedir, "BLprofiles_$(label).png")
    savefig(p, outfile)
    @info "Saved: $outfile"
    return p
end
