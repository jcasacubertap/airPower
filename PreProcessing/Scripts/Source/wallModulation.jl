# ──────────────────────────────────────────────────────────────────────
#  Wall modulation — pure shape functions (no backend / Docker deps)
#  Included by backend.jl and by generateGrid.jl.
# ──────────────────────────────────────────────────────────────────────

"""
    _smoothstep(t, n)

Generalized smoothstep: S(0)=0, S(1)=1, S'(0)=S'(1)=0 for n≥2.
"""
_smoothstep(t, n) = t^n / (t^n + (1 - t)^n)

"""
    _esn_geometry(wm) → (xStart, xEnd)

Compute the truncated base extent of the ESN bump from its parameters.
Returns (xStart, xEnd) in the same coordinates as `wm.xCenter`.
"""
function _esn_geometry(wm)
    # Unscaled ESN on wide range
    x_esn = range(-8, 8, length=10000)
    sigma_left  = 1.0 + wm.epsilon
    sigma_right = 1.0 - wm.epsilon

    y_esn = [xi < 0 ? exp(-xi^2 / (2 * sigma_left^2)) :
                       exp(-xi^2 / (2 * sigma_right^2))
             for xi in x_esn]

    # Truncate tails at yTol / A (normalized)
    tol_norm = wm.yTol / abs(wm.A)
    idx1 = findlast(i -> y_esn[i] < tol_norm && x_esn[i] < 0, 1:length(x_esn))
    idx2 = findfirst(i -> y_esn[i] < tol_norm && x_esn[i] > 0, 1:length(x_esn))
    if idx1 === nothing; idx1 = 1; end
    if idx2 === nothing; idx2 = length(x_esn); end

    b_esn = x_esn[idx2] - x_esn[idx1]  # unscaled truncated base width

    # Scale to target base width
    b_target = wm.R * abs(wm.A)
    scale = b_target / b_esn

    xStart = wm.xCenter + x_esn[idx1] * scale
    xEnd   = wm.xCenter + x_esn[idx2] * scale
    return (xStart, xEnd)
end

"""
    _esn_geometry_full(wm) → (xStart_full, xEnd_full)

Like `_esn_geometry` but includes the 5% blend zones on each side.
Used for polyLine extent and vertex checks.
"""
function _esn_geometry_full(wm)
    xStart, xEnd = _esn_geometry(wm)
    blend_w_phys = 0.05 * (xEnd - xStart)
    return (xStart - blend_w_phys, xEnd + blend_w_phys)
end

"""
    wall_bump_sigmoidal(x; A, xStart, xPeak, xEnd, p, q)

Sigmoidal bump: piecewise smoothstep, C2 for p,q ≥ 3, compact support.
"""
function wall_bump_sigmoidal(x; A, xStart, xPeak, xEnd, p, q)
    (x <= xStart || x >= xEnd) && return 0.0
    if x <= xPeak
        t = (x - xStart) / (xPeak - xStart)
        return A * _smoothstep(t, p)
    else
        s = (x - xPeak) / (xEnd - xPeak)
        return A * (1.0 - _smoothstep(s, q))
    end
end

"""
    wall_bump_esn(x; A, xCenter, epsilon, R, yTol)

Epsilon-Skewed Normal bump: split Gaussian with skewness.
Inside the yTol truncation boundaries: exact ESN shape.
Outside: C2-smooth blend from yTol to 0 over 5% of the bump width.
"""
function wall_bump_esn(x; A, xCenter, epsilon, R, yTol)
    sigma_left  = 1.0 + epsilon
    sigma_right = 1.0 - epsilon

    tol_norm = yTol / abs(A)
    tol_norm >= 1.0 && return 0.0
    x_trunc_left  = -sigma_left  * sqrt(-2.0 * log(tol_norm))
    x_trunc_right =  sigma_right * sqrt(-2.0 * log(tol_norm))
    b_esn = x_trunc_right - x_trunc_left

    b_target = R * abs(A)
    scale = b_target / b_esn

    blend_w = 0.05 * b_esn

    xi = (x - xCenter) / scale

    if xi > x_trunc_left && xi < x_trunc_right
        sigma = xi < 0 ? sigma_left : sigma_right
        return A * exp(-xi^2 / (2.0 * sigma^2))
    elseif xi >= x_trunc_left - blend_w && xi <= x_trunc_left
        sigma = sigma_left
        gauss = A * exp(-xi^2 / (2.0 * sigma^2))
        t = (xi - (x_trunc_left - blend_w)) / blend_w
        return gauss * _smoothstep(t, 3)
    elseif xi >= x_trunc_right && xi <= x_trunc_right + blend_w
        sigma = sigma_right
        gauss = A * exp(-xi^2 / (2.0 * sigma^2))
        t = (xi - x_trunc_right) / blend_w
        return gauss * (1.0 - _smoothstep(t, 3))
    else
        return 0.0
    end
end

"""
    wall_bump(x, wm)

Dispatch to the appropriate bump function based on `wm.shape`.
"""
function wall_bump(x, wm)
    if wm.shape == :sigmoidal
        return wall_bump_sigmoidal(x; A=wm.A, xStart=wm.xStart, xPeak=wm.xPeak,
                                   xEnd=wm.xEnd, p=wm.p, q=wm.q)
    elseif wm.shape == :esn
        return wall_bump_esn(x; A=wm.A, xCenter=wm.xCenter, epsilon=wm.epsilon,
                              R=wm.R, yTol=wm.yTol)
    else
        error("Unknown wall modulation shape: $(wm.shape)")
    end
end
