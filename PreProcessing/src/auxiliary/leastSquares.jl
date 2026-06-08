#
# Polynomial least-squares fitting (monomial and logarithmic bases)
# Translated from MATLAB: leastSquaresUnconstrainedFun.m, leastSquaresConstrainedFun.m
# Jordi Casacuberta, 2019 / Julia port 2026
#

"""
    build_basis_matrix(x, N, fit_law) → Matrix

Build the design matrix for polynomial fitting of order `N`.
- `fit_law = :monomial`:    columns are [1, x, x², ..., x^N]
- `fit_law = :logarithmic`: columns are [1, ln(x), ln(x)², ..., ln(x)^N]
"""
function build_basis_matrix(x::AbstractVector, N::Int, fit_law::Symbol)
    m = length(x)
    X = zeros(m, N + 1)
    basis = fit_law == :logarithmic ? log.(x) : Float64.(x)
    for j in 0:N
        X[:, j+1] .= basis .^ j
    end
    return X
end

"""
    eval_polynomial(x, coeffs, N, fit_law) → Vector

Evaluate polynomial with coefficients in descending order (highest power first).
"""
function eval_polynomial(x::AbstractVector, coeffs::AbstractVector, N::Int, fit_law::Symbol)
    basis = fit_law == :logarithmic ? log.(x) : Float64.(x)
    result = zeros(length(x))
    for i in 1:N+1
        result .+= coeffs[i] .* basis .^ (N - (i - 1))
    end
    return result
end

"""
    eval_polynomial_deriv(x_val, coeffs, N, fit_law) → Float64

First derivative of the polynomial at a single point, computed analytically.
"""
function eval_polynomial_deriv(x_val::Real, coeffs::AbstractVector, N::Int, fit_law::Symbol)
    result = 0.0
    k = 1
    if fit_law == :monomial
        for i in N:-1:1
            result += i * coeffs[k] * x_val^(i - 1)
            k += 1
        end
    elseif fit_law == :logarithmic
        for i in N:-1:1
            result += i * coeffs[k] * log(x_val)^(i - 1) / x_val
            k += 1
        end
    end
    return result
end

"""
    least_squares_unconstrained(xe, ye, N, fit_law) → coeffs

Polynomial least-squares fit (unconstrained). Returns coefficients
in descending order (highest power first).
"""
function least_squares_unconstrained(xe::AbstractVector, ye::AbstractVector,
                                     N::Int, fit_law::Symbol)
    X = build_basis_matrix(xe, N, fit_law)
    coeffs = X \ ye
    return reverse(coeffs)
end

"""
    least_squares_constrained(xe, ye, N, xt, fit_law) → coeffs

Polynomial least-squares fit with zero 1st and 2nd derivative at `xt`.
Uses Lagrange multipliers. Returns coefficients in descending order.

Source: Constrained Linear Least Squares, Henri P. Gavin, 2015, Duke Un.
"""
function least_squares_constrained(xe::AbstractVector, ye::AbstractVector,
                                    N::Int, xt::Real, fit_law::Symbol)
    X = build_basis_matrix(xe, N, fit_law)

    # Constraint matrix: A * coeffs = 0 (zero 1st and 2nd derivative at xt)
    A = zeros(2, N + 1)

    if fit_law == :monomial
        for i in 2:N+1
            A[1, i] = (i - 1) * xt^(i - 2)            # 1st derivative
        end
        for i in 3:N+1
            A[2, i] = (i - 1) * (i - 2) * xt^(i - 3)  # 2nd derivative
        end
    elseif fit_law == :logarithmic
        for i in 2:N+1
            A[1, i] = (i - 1) * log(xt)^(i - 2) / xt
        end
        A[2, 2] = -1 / xt^2
        for i in 3:N+1
            A[2, i] = -(i - 1) * (log(xt) - (i - 2)) * log(xt)^(i - 3) / xt^2
        end
    end

    b = [0.0, 0.0]

    # Augmented Lagrangian system
    M   = [2*X'*X  A';  A  zeros(2, 2)]
    rhs = [2*X'*ye; b]
    sol = M \ rhs

    coeffs = sol[1:N+1]
    return reverse(coeffs)
end
