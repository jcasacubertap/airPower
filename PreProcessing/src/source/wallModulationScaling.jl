#
# wallModulationScaling — resolve a single wall-modulation bump from any two of
# the four parameters {A, xCenter, Re_k, A/δ*}.
#
# The four are tied by two relations evaluated on the *undisturbed* baseline
# boundary layer (here supplied by the third-party IBL solver):
#
#     Re_k      = u_b(x, y=A) · A / ν          (roughness Reynolds number)
#     A/δ*      = A / δ*_b(x)                   (height over displacement thickness)
#
# with x ≡ xCenter the bump station. Both depend on (A, x) only, so the system
# has two degrees of freedom: prescribe any two of the four and the other two
# are determined. Each combination collapses to at most a 1-D, monotone
# root-find, so a bracketed bisection is robust and unique.
#
# The baseline BL enters only through `BaselineBL` (a velocity-profile evaluator
# u_b(x,y) and a displacement-thickness evaluator δ*_b(x)), so this module is
# independent of how the baseline is produced and is unit-testable on an
# analytic profile without running the solver.

"""
    BaselineBL(ub, dstar, nu, xlo, xhi, Alo, Ahi)

Undisturbed baseline boundary layer used to evaluate the bump-scaling relations.

  - `ub(x, y)`   → streamwise (chordwise) velocity at station `x`, height `y` [m/s]
  - `dstar(x)`   → displacement thickness δ*(x) [m]
  - `nu`         → kinematic viscosity [m²/s]
  - `[xlo, xhi]` → bracket for root-finds over the bump station x [m]
  - `[Alo, Ahi]` → bracket for root-finds over the bump height A [m]
"""
struct BaselineBL
    ub::Function
    dstar::Function
    nu::Float64
    xlo::Float64
    xhi::Float64
    Alo::Float64
    Ahi::Float64
end

_rek(bl::BaselineBL, A, x) = bl.ub(x, A) * A / bl.nu

# Bracketed bisection for a monotone f on [lo, hi]; clear error if the target
# lies outside the achievable range (no sign change).
function _bisect(f, lo, hi; tol = 1e-12, maxit = 200)
    flo = f(lo); fhi = f(hi)
    flo == 0 && return lo
    fhi == 0 && return hi
    sign(flo) == sign(fhi) && error(
        "wallModulationScaling: target not bracketed in [$(lo), $(hi)] " *
        "(f(lo)=$(round(flo, sigdigits=4)), f(hi)=$(round(fhi, sigdigits=4))) — " *
        "the requested value is outside the achievable range of the baseline BL.")
    mid = 0.5 * (lo + hi)
    for _ in 1:maxit
        mid = 0.5 * (lo + hi)
        fm  = f(mid)
        (fm == 0 || (hi - lo) < tol * max(1.0, abs(mid))) && return mid
        if sign(fm) == sign(flo)
            lo = mid; flo = fm
        else
            hi = mid; fhi = fm
        end
    end
    return mid
end

"""
    resolve_bump(bl; A=nothing, xCenter=nothing, Rek=nothing, AoverDstar=nothing)

Given exactly two of the four bump parameters, solve for the rest on the
baseline BL `bl`. Returns `(A, xCenter, Rek, AoverDstar)` fully populated.

Supported pairs (all 6 of `choose(4,2)`):

| given                | solved        | how                                   |
|----------------------|---------------|---------------------------------------|
| `{A, xCenter}`       | Rek, A/δ*     | forward (direct)                      |
| `{xCenter, AoverDstar}` | A, Rek     | A = (A/δ*)·δ*(x)                       |
| `{A, Rek}`           | xCenter, A/δ* | root-find x: u_b(x,A)·A/ν = Rek        |
| `{A, AoverDstar}`    | xCenter, Rek  | root-find x: δ*(x) = A/(A/δ*)          |
| `{xCenter, Rek}`     | A, A/δ*       | root-find A: u_b(x,A)·A/ν = Rek        |
| `{Rek, AoverDstar}`  | A, xCenter    | A=(A/δ*)·δ*(x), root-find x in Rek     |
"""
function resolve_bump(bl::BaselineBL; A = nothing, xCenter = nothing,
                      Rek = nothing, AoverDstar = nothing)
    spec = (; A, xCenter, Rek, AoverDstar)
    gk   = Tuple(k for k in keys(spec) if spec[k] !== nothing)
    length(gk) == 2 || error(
        "resolve_bump: specify exactly two of {A, xCenter, Rek, AoverDstar}; got $(gk).")

    s = Set(gk)
    a = A; x = xCenter; r = AoverDstar

    if s == Set((:A, :xCenter))
        # forward — nothing to solve
    elseif s == Set((:xCenter, :AoverDstar))
        a = r * bl.dstar(x)
    elseif s == Set((:A, :Rek))
        x = _bisect(xx -> _rek(bl, a, xx) - Rek, bl.xlo, bl.xhi)
    elseif s == Set((:A, :AoverDstar))
        x = _bisect(xx -> bl.dstar(xx) - a / r, bl.xlo, bl.xhi)
    elseif s == Set((:xCenter, :Rek))
        a = _bisect(aa -> _rek(bl, aa, x) - Rek, bl.Alo, bl.Ahi)
    elseif s == Set((:Rek, :AoverDstar))
        x = _bisect(xx -> _rek(bl, r * bl.dstar(xx), xx) - Rek, bl.xlo, bl.xhi)
        a = r * bl.dstar(x)
    else
        error("resolve_bump: unsupported combination $(gk).")
    end

    return (A = a, xCenter = x, Rek = _rek(bl, a, x), AoverDstar = a / bl.dstar(x))
end

# ── Baseline BL from the third-party IBL solution ─────────────────────────────
# Small self-contained interpolators (no external deps) so the adapter is
# testable on a synthetic IBL-shaped solution without MATLAB.

_trapz(x, y) = sum((x[i+1] - x[i]) * (y[i+1] + y[i]) / 2 for i in 1:length(x)-1)

function _interp1(xs, ys, x)
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs, x)
    t = (x - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] + t * (ys[k+1] - ys[k])
end

# bilinear on an ascending grid; U indexed [ix, iy]
function _bilinear(xs, ys, U, x, y)
    nx = length(xs); nyv = length(ys)
    kx = clamp(searchsortedlast(xs, x), 1, nx - 1)
    ky = clamp(searchsortedlast(ys, y), 1, nyv - 1)
    tx = clamp((x - xs[kx]) / (xs[kx+1] - xs[kx]), 0.0, 1.0)
    ty = clamp((y - ys[ky]) / (ys[ky+1] - ys[ky]), 0.0, 1.0)
    return (1-tx)*(1-ty)*U[kx,ky]   + tx*(1-ty)*U[kx+1,ky] +
           (1-tx)*ty  *U[kx,ky+1]   + tx*ty    *U[kx+1,ky+1]
end

"""
    baseline_from_ibl(sol; nu, Amax=nothing, xtrim=0.0) → BaselineBL

Build a `BaselineBL` from an IBL solution `sol = (Xgrid, Y, u, …)` (the
NamedTuple returned by `run_ibl`/`run_ibl_solver`). `u` is the chordwise
velocity (ny × nx); `Y` the wall-normal coordinate; `Xgrid` the inlet-referenced
stations (same frame as the DFP `xCenter`). δ*(x) is `∫(1−u/U_e)dy` with the
edge velocity `U_e` taken as the profile's outermost value. `xtrim` shrinks the
x bracket by that fraction at each end; `Amax` caps the A bracket (defaults to
the domain top).
"""
function baseline_from_ibl(sol; nu::Real, Amax::Union{Real,Nothing}=nothing,
                           xtrim::Real=0.0)
    Xg = collect(Float64.(sol.Xgrid))
    Yv = collect(Float64.(sol.Y))
    U  = Matrix{Float64}(sol.u)                 # ny × nx (rows = y, cols = x)

    px = sortperm(Xg); Xg = Xg[px]; U = U[:, px]
    py = sortperm(Yv); Yv = Yv[py]; U = U[py, :]   # rows now ascending in y (wall→edge)
    Uxy = permutedims(U)                          # nx × ny, indexed [ix, iy]
    nx  = length(Xg)

    ds = Vector{Float64}(undef, nx)
    for k in 1:nx
        prof = Uxy[k, :]                          # ascending y, wall→edge
        ds[k] = _trapz(Yv, 1 .- prof ./ prof[end])
    end

    ub(x, y) = _bilinear(Xg, Yv, Uxy, x, y)
    dstar(x) = _interp1(Xg, ds, x)
    span = Xg[end] - Xg[1]
    return BaselineBL(ub, dstar, Float64(nu),
                      Xg[1] + xtrim*span, Xg[end] - xtrim*span,
                      1e-7, Amax === nothing ? Yv[end] : Float64(Amax))
end

"""
    ibl_baseline_bl(; kwargs...) → BaselineBL

Run the third-party IBL solver for the current DFP inputs and wrap it as a
`BaselineBL`. Requires `run_ibl_solver` and `inp` (loaded via run.jl).
"""
function ibl_baseline_bl(; kwargs...)
    sol = run_ibl_solver()
    sol === nothing && error("ibl_baseline_bl: IBL baseline unavailable (MATLAB/solver failed).")
    return baseline_from_ibl(sol; nu = inp.DFP.freeStreamViscosity, kwargs...)
end

# ── DFP wiring: resolve the bump geometry from the inputs spec (cached) ───────

const _DFP_BUMP_CACHE = Ref{Union{Nothing,NamedTuple}}(nothing)

"Clear the cached DFP bump resolution (e.g. after changing inputs in a REPL)."
reset_dfp_bump!() = (_DFP_BUMP_CACHE[] = nothing)

"""
    dfp_bump_geometry() → (A, xCenter, Rek, AoverDstar)

Resolve the DFP bump from `inp.DFP.wallModulation`. Specify exactly two of
{A, xCenter, Rek, AoverDstar}: a directly-given {A, xCenter} skips the solver;
any Rek/AoverDstar target runs the IBL baseline once (cached) and `resolve_bump`
solves the rest. The resolved `A` overrides the shared `inp.wallModulation.A`
for the DFP. (ESN/`xCenter` only for now; sigmoidal scaling is future work.)
"""
function dfp_bump_geometry()
    _DFP_BUMP_CACHE[] !== nothing && return _DFP_BUMP_CACHE[]
    s   = inp.DFP.wallModulation
    A   = hasproperty(s, :A)          ? s.A          : nothing
    xc  = hasproperty(s, :xCenter)    ? s.xCenter    : nothing
    Rek = hasproperty(s, :Rek)        ? s.Rek        : nothing
    aod = hasproperty(s, :AoverDstar) ? s.AoverDstar : nothing

    res = if Rek !== nothing || aod !== nothing
        count(!isnothing, (A, xc, Rek, aod)) == 2 || error(
            "DFP.wallModulation: with a scaling target, give exactly two of " *
            "{A, xCenter, Rek, AoverDstar}.")
        @info "Resolving DFP bump from scaling targets via the IBL baseline…"
        r = resolve_bump(ibl_baseline_bl(); A = A, xCenter = xc, Rek = Rek, AoverDstar = aod)
        @info "DFP bump resolved" A = r.A xCenter = r.xCenter Rek = round(r.Rek, digits = 2) AoverDstar = round(r.AoverDstar, digits = 4)
        r
    else
        xc === nothing && error("DFP.wallModulation: direct mode needs xCenter (with A, or shared wallModulation.A).")
        Aeff = A !== nothing ? A : inp.wallModulation.A
        (A = Aeff, xCenter = xc, Rek = nothing, AoverDstar = nothing)
    end
    _DFP_BUMP_CACHE[] = res
    return res
end

"""
    dfp_wm() → NamedTuple

DFP wall-modulation parameters (shared shape ∪ DFP position) with `A`/`xCenter`
set to the resolved bump geometry. Use everywhere the DFP needs the bump
(`write_flat_plate_input_param`, mesh, viz) so direct and scaling specs behave
identically downstream.
"""
function dfp_wm()
    wm = merge(inp.wallModulation, inp.DFP.wallModulation)
    if wm.enabled
        g  = dfp_bump_geometry()
        wm = merge(wm, (A = g.A, xCenter = g.xCenter))
    end
    return wm
end
