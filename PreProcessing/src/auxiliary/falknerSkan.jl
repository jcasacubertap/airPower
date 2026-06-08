#
# Falkner-Skan-Cooke similarity solver (shooting + RK4)
# Translated from MATLAB: FalknerSkanSolverFun.m, RK4.m
# Jordi Casacuberta, 2019 / Julia port 2026
#

"""
    rk4_falkner_skan!(u1, u2, u3, u4, u5, h, beta)

Fourth-order Runge-Kutta integrator for the Falkner-Skan-Cooke system:

    f'''  = -f f'' - beta (1 - f'^2)
    g''   = -f g'

where u1=f, u2=f', u3=f'', u4=g, u5=g'.
Marches in-place from index 1 to length(u1).
"""
function rk4_falkner_skan!(u1, u2, u3, u4, u5, h::Real, beta::Real)
    N = length(u1)
    for i in 1:N-1
        # Stage 1
        k1 = h * u2[i]
        l1 = h * u3[i]
        r1 = h * (-u1[i] * u3[i] - beta * (1 - u2[i]^2))
        n1 = h * u5[i]
        t1 = h * (-u1[i] * u5[i])

        # Stage 2
        k2 = h * (u2[i] + l1 / 2)
        l2 = h * (u3[i] + r1 / 2)
        r2 = h * (-(u1[i] + k1 / 2) * (u3[i] + r1 / 2) - beta * (1 - (u2[i] + l1 / 2)^2))
        n2 = h * (u5[i] + t1 / 2)
        t2 = h * (-(u1[i] + k1 / 2) * (u5[i] + t1 / 2))

        # Stage 3
        k3 = h * (u2[i] + l2 / 2)
        l3 = h * (u3[i] + r2 / 2)
        r3 = h * (-(u1[i] + k2 / 2) * (u3[i] + r2 / 2) - beta * (1 - (u2[i] + l2 / 2)^2))
        n3 = h * (u5[i] + t2 / 2)
        t3 = h * (-(u1[i] + k2 / 2) * (u5[i] + t2 / 2))

        # Stage 4
        k4 = h * (u2[i] + l3)
        l4 = h * (u3[i] + r3)
        r4 = h * (-(u1[i] + k3) * (u3[i] + r3) - beta * (1 - (u2[i] + l3)^2))
        n4 = h * (u5[i] + t3)
        t4 = h * (-(u1[i] + k3) * (u5[i] + t3))

        # Update
        u1[i+1] = u1[i] + (k1 + 2k2 + 2k3 + k4) / 6
        u2[i+1] = u2[i] + (l1 + 2l2 + 2l3 + l4) / 6
        u3[i+1] = u3[i] + (r1 + 2r2 + 2r3 + r4) / 6
        u4[i+1] = u4[i] + (n1 + 2n2 + 2n3 + n4) / 6
        u5[i+1] = u5[i] + (t1 + 2t2 + 2t3 + t4) / 6
    end
    return u1, u2, u3, u4, u5
end

"""
    falkner_skan_solve(beta_FSK) → (u1, u2, u3, u4, u5, eta, m)

Solve the Falkner-Skan-Cooke equation for a given pressure-gradient
parameter beta. Uses a shooting method with secant iteration.

Returns: f, f', f'', g, g', eta, m (Hartree parameter).
"""
function falkner_skan_solve(beta_FSK::Real)
    m = beta_FSK / (2 - beta_FSK)
    maxetaglobal = 200.0

    etaf_init = (-0.05 ≤ beta_FSK ≤ 1.5) ? 3.8 : 2.8

    h       = 1e-2
    tol     = 1e-12
    itermax = 100

    # --- Secant shooting on a given domain ----------------------------------
    function shoot(etaf, alpha1, alpha2, theta1, theta2)
        N = round(Int, etaf / h)
        u1 = zeros(N); u2 = zeros(N); u3 = zeros(N)
        u4 = zeros(N); u5 = zeros(N)
        alpha  = zeros(itermax); solu2 = zeros(itermax)
        theta  = zeros(itermax); solu4 = zeros(itermax)

        # Two conditioning runs
        for run in 1:2
            fill!(u1, 0.0); fill!(u2, 0.0); fill!(u3, 0.0)
            fill!(u4, 0.0); fill!(u5, 0.0)
            alpha[run] = run == 1 ? alpha1 : alpha2
            theta[run] = run == 1 ? theta1 : theta2
            u3[1] = alpha[run]
            u5[1] = theta[run]
            rk4_falkner_skan!(u1, u2, u3, u4, u5, h, beta_FSK)
            solu2[run] = u2[end]
            solu4[run] = u4[end]
        end

        # Secant iteration
        for j in 3:itermax
            fill!(u1, 0.0); fill!(u2, 0.0); fill!(u3, 0.0)
            fill!(u4, 0.0); fill!(u5, 0.0)

            slope_a   = (alpha[j-1] - alpha[j-2]) / (solu2[j-1] - solu2[j-2])
            alpha[j]  = alpha[j-1] - slope_a * (solu2[j-1] - 1)
            u3[1]     = alpha[j]

            slope_t   = (theta[j-1] - theta[j-2]) / (solu4[j-1] - solu4[j-2])
            theta[j]  = theta[j-1] - slope_t * (solu4[j-1] - 1)
            u5[1]     = theta[j]

            rk4_falkner_skan!(u1, u2, u3, u4, u5, h, beta_FSK)
            solu2[j] = u2[end]
            solu4[j] = u4[end]

            if abs(solu2[j] - 1) < tol && abs(solu4[j] - 1) < tol
                break
            end
        end
        return u3[1], u5[1]
    end

    # Phase 1: small domain — find initial guesses for f''(0) and g'(0)
    val3, val5 = shoot(etaf_init, 1.3, 1.35, 1.1, 1.2)

    # Phase 2: medium domain (etaf = 10) — refine
    val33, val55 = shoot(10.0, val3, val3 + 1e-4, val5, val5 + 1e-4)

    # Phase 3: final integration on the full domain
    N   = round(Int, maxetaglobal / h)
    eta = collect(range(0.0, maxetaglobal, length=N))
    u1  = zeros(N); u2 = zeros(N); u3 = zeros(N)
    u4  = zeros(N); u5 = zeros(N)
    u3[1] = val33
    u5[1] = val55
    rk4_falkner_skan!(u1, u2, u3, u4, u5, h, beta_FSK)

    return u1, u2, u3, u4, u5, eta, m
end
