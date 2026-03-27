#
# External flow data → scaling parameters for the flat-plate computation
# Jordi Casacuberta, 2019 / Julia port 2026
#
# Steps:
#   1. Crop raw external data (x, ue) from a chord % to monotonicity loss
#   2. Fit a polynomial P(x) on the cropped data (monomial or log basis)
#   3. Derive scaling: uinf = P(x_crop[1]), due/dx, beta_FSK, delta0
#

using MAT, Interpolations

"""
    external_to_scaling(; x, S, expUe, chord, percent_crop, Nfit, fit_law,
                          nuinf, savedir) → NamedTuple

Process an external edge-velocity distribution and return scaling parameters
for the flat-plate computation.

# Arguments
- `x`:            streamwise coordinate (from airfoil geometry) [m]
- `S`:            arc-length coordinate [m]
- `expUe`:        experimental/external edge velocity [m/s]
- `chord`:        airfoil chord [m]
- `percent_crop`: crop start as % of chord
- `Nfit`:         polynomial order
- `fit_law`:      `:monomial` or `:logarithmic`
- `nuinf`:        free-stream kinematic viscosity [m^2/s]
- `savedir`:      directory for diagnostic plots (nothing to skip)
"""
function external_to_scaling(; x::AbstractVector, S::AbstractVector,
                               expUe::AbstractVector, chord::Real,
                               percent_crop::Real, Nfit::Int,
                               fit_law::Symbol, nuinf::Real,
                               savedir = nothing)

    # ── 1. Crop raw data ────────────────────────────────────────────────────
    # Start: first point where x >= percent_crop/100 * chord
    st = findfirst(xi -> xi >= percent_crop / 100 * chord, x)

    # End: streamwise loss of monotonicity in expUe
    up_lim = length(expUe)
    for i in 2:length(expUe)
        if expUe[i] < expUe[i-1]
            up_lim = i - 1
            break
        end
    end

    expUe_crop = expUe[st:up_lim]
    S_crop     = S[st:up_lim]
    xt         = S[up_lim]

    # ── 2. Polynomial fit P(x) on cropped data ─────────────────────────────
    ue_pol_coeff = least_squares_unconstrained(S_crop, expUe_crop, Nfit, fit_law)
    ue_pol       = eval_polynomial(S_crop, ue_pol_coeff, Nfit, fit_law)

    # Characteristic velocity: edge velocity at the inlet
    uinf = ue_pol[1]

    # ── 3. Falkner-Skan-Cooke scaling ──────────────────────────────────────
    # Analytical due/dx at the inlet from the polynomial
    duedx_1 = eval_polynomial_deriv(S_crop[1], ue_pol_coeff, Nfit, fit_law)

    # Finite-difference check (2nd-order forward)
    duedx_1_check = (-3ue_pol[1] + 4ue_pol[2] - ue_pol[3]) /
                    (2 * (S_crop[2] - S_crop[1]))

    # Falkner-Skan and Hartree parameters
    m_hartree = duedx_1 * S_crop[1] / ue_pol[1]
    beta_FSK  = 2 * m_hartree / (1 + m_hartree)

    # Solve Falkner-Skan-Cooke similarity equation
    _, u2_fsk, _, _, _, eta_fsk, _ = falkner_skan_solve(beta_FSK)

    # Rescale to physical coordinates (flip: free-stream → wall)
    y_re = reverse(eta_fsk) ./ sqrt((m_hartree + 1) / 2 *
           (ue_pol[1] / (S_crop[1] * nuinf)))
    u_re = reverse(ue_pol[1] .* u2_fsk)

    # Find delta_99: where u first drops below 0.99 * ue
    u99 = 0.99 * ue_pol[1]
    k   = findfirst(u -> u <= u99, u_re)

    # Interpolation around the crossing for delta0
    idx_range = max(1, k - 3):min(length(u_re), k + 3)
    u_local   = u_re[idx_range]
    y_local   = y_re[idx_range]
    perm      = sortperm(u_local)
    itp       = linear_interpolation(u_local[perm], y_local[perm])
    delta0    = itp(u99)

    # ── 4. Diagnostic plots ─────────────────────────────────────────────────
    if savedir !== nothing
        _plot_external_scaling(; S, expUe, S_crop, expUe_crop, ue_pol, xt,
                                 u_re, y_re, delta0, uinf,
                                 x_geom=x, st, percent_crop, savedir)
    end

    # ── Return ──────────────────────────────────────────────────────────────
    return (
        delta0        = delta0,
        uinf          = uinf,
        beta_FSK      = beta_FSK,
        duedx_1       = duedx_1,
        duedx_1_check = duedx_1_check,
        m_hartree     = m_hartree,
        S_crop        = S_crop,
        expUe_crop    = expUe_crop,
        xt            = xt,
        ue_pol_coeff  = ue_pol_coeff,
        ue_pol        = ue_pol,
    )
end

# ── Diagnostic plots ────────────────────────────────────────────────────────
function _plot_external_scaling(; S, expUe, S_crop, expUe_crop, ue_pol, xt,
                                  u_re, y_re, delta0, uinf,
                                  x_geom, st, percent_crop, savedir)
    mkpath(savedir)

    # (a) Dimensional edge velocity: raw, cropped, polynomial fit
    p1 = plot(S, expUe, color=:black, lw=1.2, label="raw data")
    plot!(p1, S_crop, expUe_crop, color=:green, lw=1.2,
          markershape=:circle, ms=2, label="cropped raw data")
    plot!(p1, S_crop, ue_pol, color=:red, lw=1.2,
          markershape=:circle, ms=2, label="polynomial fit")
    vline!(p1, [xt], color=:black, ls=:dash, label=nothing)
    xlabel!(p1, L"S_{\mathrm{exp}}\;(\mathrm{m})")
    ylabel!(p1, L"u_e\;(\mathrm{m/s})")
    ylims!(p1, (0, maximum(expUe) * 1.2))

    # (b) Inflow BL profile from Falkner-Skan
    p2 = plot(u_re, y_re, color=:black, lw=1.2, label=nothing)
    hline!(p2, [delta0], color=:black, ls=:dash, label=nothing)
    xlabel!(p2, L"u_{\mathrm{inflow}}\;(\mathrm{m/s})")
    ylabel!(p2, L"y\;(\mathrm{m})")
    ylims!(p2, (0, 3delta0))

    # (c) Geometrical relation x → S
    p3 = plot(x_geom, S, color=:black, lw=1.2,
              markershape=:circle, ms=2, label="geometrical relation")
    vline!(p3, [x_geom[st]], color=:magenta, ls=:dash,
           label="$(percent_crop)% chord")
    xlabel!(p3, L"x_{\mathrm{exp}}\;(\mathrm{m})")
    ylabel!(p3, L"S_{\mathrm{exp}}\;(\mathrm{m})")
    ylims!(p3, (0, maximum(S)))

    fig = plot(p1, p2, p3, layout=(2, 2), size=(900, 700),
               plot_title="Reference external flow properties")

    savefig(fig, joinpath(savedir, "externalToScaling.png"))
    @info "Saved diagnostic plot" path=joinpath(savedir, "externalToScaling.png")
    return fig
end
