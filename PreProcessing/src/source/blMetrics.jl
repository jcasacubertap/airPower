#!/usr/bin/env julia
#
# blMetrics — boundary-layer integral metrics post-processor.
#
# For each suction-side wall face exported by the writeMidPlane function
# object, compute (d99, dstar, Theta) and the freestream speed U_e at that
# station, and write postProcessing/BLQuantities.csv (columns x, s, xi, Ue,
# d99, dstar, Theta — x is the global, post-AoA-rotation wall-face x in
# metres (same convention as postProcessing/wallQuantities.csv), s is the
# arclength along the upper surface measured from xi=0, and xi is the chord
# fraction (x/c) used to plot vs the airfoil-chord coordinate).
# Each mid-plane cell is Voronoi-assigned to its nearest wall face to build
# the wall-normal column.
#
# The freestream method is selected by `blMetrics.method` in
# constant/inputParam. Currently supported:
#   vorticityIntegralTrapezoidal — U_e = −∫ ω_z dn along the wall normal,
#                       evaluated as a trapezoidal sum over interior mid-plane
#                       cells (sorted ascending in wall-normal distance n).
#                       Near-wall handling: the wall-extrapolation row written
#                       by writeMidPlane (a single (x_w, y_w, z, 0, 0, 0, p,
#                       omz) tuple at n ≈ 0, with ω_z linearly extrapolated
#                       from the two innermost interior cells) is *excluded*
#                       from the integral, because its extrapolated ω_z can be
#                       unphysical near the leading edge. As a consequence the
#                       integral starts at n_int[1] (the first interior cell,
#                       typically ~10–100 μm above the wall) and the
#                       contribution from n=0 to n_int[1] — where |ω_z| peaks
#                       — is missed. This biases U_e slightly low, especially
#                       in regions of strong wall-normal vorticity (e.g., near
#                       the LE). The other methods (maxProfile, fixedHeight)
#                       are immune to this near-wall slice.
#
# Invocation:
#   julia blMetrics.jl <caseDir>          # caseDir defaults to "."
#
# Jordi Casacuberta, 2026

using DelimitedFiles
using Printf
using Statistics

# ---------------------------------------------------------------------------
# Tiny OpenFOAM dictionary parser (only what blMetrics needs)
# ---------------------------------------------------------------------------

"""
    read_subdict(path, name) -> Dict{String,Any}

Parse a top-level sub-dictionary `<name> { ... }` from an OpenFOAM dict file.
Returns scalars as Float64 (or String if not numeric) and `N(v1 v2 ...)` lists
as Vector{Float64}. Comments (`//` and `/* */`) are stripped.
"""
function read_subdict(path::AbstractString, name::AbstractString)
    txt = read(path, String)
    txt = replace(txt, r"/\*[\s\S]*?\*/" => " ")
    txt = replace(txt, r"//[^\n]*" => "")

    m = match(Regex("\\b" * name * "\\s*\\{", "s"), txt)
    m === nothing && error("Sub-dictionary '$(name)' not found in $(path)")

    # Find the matching closing brace
    i = m.offset + length(m.match)
    depth = 1; j = i
    while j <= lastindex(txt) && depth > 0
        c = txt[j]
        c == '{' && (depth += 1)
        c == '}' && (depth -= 1; depth == 0 && break)
        j = nextind(txt, j)
    end
    body = txt[i:prevind(txt, j)]

    out = Dict{String,Any}()
    pos = firstindex(body)
    while true
        while pos <= lastindex(body) && isspace(body[pos])
            pos = nextind(body, pos)
        end
        pos > lastindex(body) && break

        # key
        ks = pos
        while pos <= lastindex(body) && !isspace(body[pos]) &&
              body[pos] != ';' && body[pos] != '{'
            pos = nextind(body, pos)
        end
        key = body[ks:prevind(body, pos)]

        while pos <= lastindex(body) && isspace(body[pos])
            pos = nextind(body, pos)
        end
        pos > lastindex(body) && break

        if body[pos] == '{'
            # Skip any nested dict (none expected at this level today)
            d = 1; pos = nextind(body, pos)
            while pos <= lastindex(body) && d > 0
                body[pos] == '{' && (d += 1)
                body[pos] == '}' && (d -= 1)
                pos = nextind(body, pos)
            end
            continue
        end

        # value up to ';'
        vs = pos
        while pos <= lastindex(body) && body[pos] != ';'
            pos = nextind(body, pos)
        end
        val = strip(body[vs:prevind(body, pos)])
        pos <= lastindex(body) && (pos = nextind(body, pos))

        list_m = match(r"^(\d+)\s*\(([^)]*)\)\s*$"s, val)
        if list_m !== nothing
            n = parse(Int, list_m.captures[1])
            toks = split(list_m.captures[2])
            length(toks) == n ||
                @warn "List size mismatch for $(key): declared $n got $(length(toks))"
            out[key] = parse.(Float64, toks)
        else
            sval = strip(val)
            f = tryparse(Float64, sval)
            out[key] = f === nothing ? String(sval) : f
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# CSV readers
# ---------------------------------------------------------------------------

"Read midPlane.csv. Returns (header::Vector{String}, data::Matrix{Float64})."
function read_midplane(path::AbstractString)
    data, hdr = readdlm(path, ',', Float64, '\n'; header=true)
    return vec(String.(hdr)), data
end

"Read wallCoordinates.csv. Returns (x::Vector, y::Vector)."
function read_wall_coords(path::AbstractString)
    data, _ = readdlm(path, ',', Float64, '\n'; header=true)
    return data[:, 1], data[:, 2]
end

# ---------------------------------------------------------------------------
# Geometry — mirrors the C++ helpers used in the writeMidPlane FO so the
# Julia post-processor and the OpenFOAM cropping logic agree on every wall
# point (same chord scaling → AoA rotation → centering).
# ---------------------------------------------------------------------------

struct Geom
    XI::Vector{Float64}      # xi/c table (LE → TE), normalised
    ETA::Vector{Float64}     # eta/c table aligned with XI (upper surface)
    chordM::Float64          # chord in metres
    alpha::Float64           # AoA in radians (signed)
    xCtrM::Float64           # rotation/translation centre x [m]
    yCtrM::Float64           # rotation/translation centre y [m]
    S::Vector{Float64}       # cumulative arclength along upper surface at XI nodes [m]
end

"Linear interpolation in a strictly increasing table xd, with clamping at ends."
@inline function interp1(xd::AbstractVector, yd::AbstractVector, xq::Real)
    n = length(xd)
    xq <= xd[1]   && return yd[1]
    xq >= xd[end] && return yd[end]
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        xd[mid] <= xq ? (lo = mid) : (hi = mid)
    end
    t = (xq - xd[lo]) / (xd[hi] - xd[lo])
    return (1 - t) * yd[lo] + t * yd[hi]
end

function build_geom(up::Dict)
    XI    = Vector{Float64}(up["airfoilXi"])
    ETA   = Vector{Float64}(up["airfoilEta"])
    chord = up["chord"]    * 1e-3       # mm → m
    alpha = up["alphaDeg"] * pi / 180   # deg → rad
    xCtr  = up["xCenter"]  * 1e-3
    yCtr  = up["yCenter"]  * 1e-3

    # Arclength along the upper surface in the pre-rotation chord-aligned
    # frame; rotation is rigid so it preserves arclengths.
    S = similar(XI); S[1] = 0.0
    @inbounds for i in 2:length(XI)
        dx = (XI[i]  - XI[i-1])  * chord
        dy = (ETA[i] - ETA[i-1]) * chord
        S[i] = S[i-1] + hypot(dx, dy)
    end
    return Geom(XI, ETA, chord, alpha, xCtr, yCtr, S)
end

eta_upper(g::Geom, xi::Real) = interp1(g.XI,  g.ETA, xi)
arclen(g::Geom,    xi::Real) = interp1(g.XI,  g.S,   xi)
xi_of_s(g::Geom,   s::Real)  = interp1(g.S,   g.XI,  s)

"Transform local (xi, eta) → global (x, y) [m]: chord-scale, rotate by α, translate."
function to_global(g::Geom, xi::Real, eta::Real)
    ca, sa = cos(g.alpha), sin(g.alpha)
    Xp =  xi  * g.chordM - 0.5 * g.chordM
    Yp =  eta * g.chordM
    Xr =  ca * Xp + sa * Yp
    Yr = -sa * Xp + ca * Yp
    return g.xCtrM + Xr, g.yCtrM + Yr
end

"Outward unit normal and unit tangent (toward increasing xi) on the upper surface at xi."
function surface_frame(g::Geom, xi::Real)
    δ = 1e-5
    xi0 = clamp(xi - δ, 0.0, 1.0)
    xi1 = clamp(xi + δ, 0.0, 1.0)
    x0, y0 = to_global(g, xi0, eta_upper(g, xi0))
    x1, y1 = to_global(g, xi1, eta_upper(g, xi1))
    tx, ty = x1 - x0, y1 - y0
    mt = hypot(tx, ty); tx /= mt; ty /= mt
    nx, ny = -ty, tx        # +90° from tangent → outward (upper surface)
    return (nx, ny, tx, ty)
end

surface_point(g::Geom, xi::Real) = to_global(g, xi, eta_upper(g, xi))

"""
    project_wall_to_xi(g, xw, yw; nref=4001) -> (xis, ss)

For each wall face at (xw[j], yw[j]), find its chord fraction xi and arclength s
by snapping to the nearest point on a refined upper-surface curve.
"""
function project_wall_to_xi(g::Geom, xw::AbstractVector, yw::AbstractVector;
                            nref::Int=4001)
    xi_ref = collect(range(0.0, 1.0; length=nref))
    XY = [to_global(g, ξ, eta_upper(g, ξ)) for ξ in xi_ref]
    S  = [arclen(g, ξ) for ξ in xi_ref]
    n  = length(xw)
    xis = Vector{Float64}(undef, n)
    ss  = Vector{Float64}(undef, n)
    @inbounds for j in 1:n
        bestk = 1; bestd = Inf
        for k in 1:nref
            d = (xw[j] - XY[k][1])^2 + (yw[j] - XY[k][2])^2
            if d < bestd
                bestd = d; bestk = k
            end
        end
        xis[j] = xi_ref[bestk]
        ss[j]  = S[bestk]
    end
    return xis, ss
end

# ---------------------------------------------------------------------------
# Cell-to-wall assignment (Voronoi): each interior/wall-extrap row of
# midPlane.csv is bound to its nearest wall face in 2D Euclidean distance.
# ---------------------------------------------------------------------------

"For each mid-plane row (xc[i], yc[i]) return the index of the nearest wall face."
function assign_cells_to_walls(xc::AbstractVector, yc::AbstractVector,
                               xw::AbstractVector, yw::AbstractVector)
    nC = length(xc); nW = length(xw)
    near = Vector{Int}(undef, nC)
    @inbounds for i in 1:nC
        bestj = 1; bestd = Inf
        xi = xc[i]; yi = yc[i]
        for j in 1:nW
            d = (xi - xw[j])^2 + (yi - yw[j])^2
            if d < bestd
                bestd = d; bestj = j
            end
        end
        near[i] = bestj
    end
    return near
end

# ---------------------------------------------------------------------------
# Per-face BL integrals
# ---------------------------------------------------------------------------

"""
    METHODS

Constants for the supported U_e (freestream-edge speed) methods.
"""
const METHODS = ("vorticityIntegralTrapezoidal", "vorticityIntegralMidpoint",
                 "maxProfile", "fixedHeight", "pressureBernoulli")

"""
    compute_ue(n_int, u_int, omz_int, method) → U_e (Float64 or NaN)

Four simple methods, all evaluated on the interior-cell column sorted
ascending in `n_int`. The wall-extrapolation row (n ≈ 0, written by the
writeMidPlane FO with u=v=w=0 and an extrapolated ω_z) has already been
filtered out by the caller, so `n_int[1]` is the first interior mid-plane
cell above the wall (≈ 10–100 μm).

  - `"vorticityIntegralTrapezoidal"` :
        U_e = − ∫_{n_int[1]}^{n_int[end]} ω_z(n) dn
        Trapezoidal sum *between cell centres*. Cell 1's weight is only
        ½·(n_int[2] − n_int[1]); the strip [0, n_int[1]] where |ω_z| peaks
        sits outside the integration interval.
  - `"vorticityIntegralMidpoint"` :
        FV-consistent midpoint rule. Each cell's centre value is multiplied
        by the cell's true wall-normal extent (face-to-face). Cell 1's
        extent extends down to the wall — it is ½·(n_int[1] + n_int[2]),
        so the wall slice is naturally included with no extrapolation.
        Interior cell weights are identical to the trapezoidal rule; the
        difference vs `vorticityIntegralTrapezoidal` is exactly the wall
        slice ω_z[1]·n_int[1] (plus a negligible top-cell adjustment).
  - `"maxProfile"`        : U_e = max(u_tan) along the column.
  - `"fixedHeight"`       : U_e = u_tan at the topmost cell (≈ exportHeight).
  - `"pressureBernoulli"` : classical Prandtl/Bernoulli.
        Assume p is constant across the BL (so the edge pressure equals the
        wall pressure: p_e = p_wall), then apply Bernoulli on a streamline
        reaching freestream conditions (U_inf, p_∞):
          U_e = sqrt( U_inf² + 2·(p_∞ − p_wall) )      (kinematic p)
        Implemented in `main()` (needs p_wall from wallQuantities.csv),
        not in `compute_ue`. Assumes p_∞ = 0 and U_inf is the chordwise
        free-stream speed read from `Upinlet.freeStreamVelocityStreamwise`.
"""
function compute_ue(n_int::Vector{Float64}, u_int::Vector{Float64},
                    omz_int::Vector{Float64}, method::AbstractString)
    K = length(n_int)
    K < 2 && return NaN
    if method == "vorticityIntegralTrapezoidal"
        # Trapezoidal integral of −ω_z from the first interior cell upwards.
        # The wall slice [0, n_int[1]] is omitted (see docstring).
        Ue = 0.0
        @inbounds for k in 1:K-1
            Ue += -0.5 * (omz_int[k] + omz_int[k+1]) * (n_int[k+1] - n_int[k])
        end
        return Ue
    elseif method == "vorticityIntegralMidpoint"
        # Per-cell midpoint rule: U_e = − Σ ω_z[k] · Δn_k where Δn_k is the
        # face-to-face extent of cell k. Face positions inferred as midpoints
        # between adjacent centres; cell 1 anchored at the wall (n=0).
        Ue = 0.0
        Ue += -omz_int[1] * 0.5 * (n_int[1] + n_int[2])             # wall → face_{1,2}
        @inbounds for k in 2:K-1
            Ue += -omz_int[k] * 0.5 * (n_int[k+1] - n_int[k-1])     # face_{k-1,k} → face_{k,k+1}
        end
        Ue += -omz_int[K] * 0.5 * (n_int[K] - n_int[K-1])           # same as trap at top (ω_z[K] ≈ 0)
        return Ue
    elseif method == "maxProfile"
        return maximum(u_int)
    elseif method == "fixedHeight"
        return u_int[end]            # topmost cell (column sorted ascending in n)
    else
        error("Unknown freestream method: $(method) — supported: $(METHODS)")
    end
end

"""
    compute_bl_integrals(n_int, u_int, Ue; max_n) → (δ99, δ*, θ)

For one wall face and a given freestream `Ue`, compute the BL thicknesses.
Anchors u_tan(0) = 0 at the wall (no slip); integrates trapezoidally up to
δ99 (first crossing of u_tan / U_e = 0.99, linearly interpolated). Returns
NaNs if the column is too short, U_e is non-positive, or δ99 is not reached.

Note: the integration is bounded by δ99 rather than extending to the top of
the column. Extending into the freestream amplifies the (1 − u/Ue) error
across a much wider tail than the BL itself when Ue is even slightly off
(1–3 % Ue error becomes ≈ 100 % δ* error in 30·δ99-thick columns).
"""
function compute_bl_integrals(n_int::Vector{Float64}, u_int::Vector{Float64},
                              Ue::Real; max_n::Real)
    K = length(n_int)
    K < 2 && return (NaN, NaN, NaN)
    (!(Ue > 0) || !isfinite(Ue)) && return (NaN, NaN, NaN)

    # Sort + truncate at max_n
    perm = sortperm(n_int)
    n_int = n_int[perm]; u_int = u_int[perm]
    if isfinite(max_n)
        cut = findlast(<=(max_n), n_int)
        cut === nothing && return (NaN, NaN, NaN)
        if cut < K
            n_int = n_int[1:cut]; u_int = u_int[1:cut]
        end
    end
    length(n_int) < 2 && return (NaN, NaN, NaN)

    # Wall point (n=0, u=0)
    n_aug = pushfirst!(copy(n_int), 0.0)
    u_aug = pushfirst!(copy(u_int), 0.0)

    # δ99 — first u/Ue = 0.99 crossing
    δ99 = NaN
    @inbounds for k in 2:length(n_aug)
        r1 = u_aug[k-1] / Ue; r2 = u_aug[k] / Ue
        if r1 < 0.99 <= r2
            t = (0.99 - r1) / (r2 - r1)
            δ99 = n_aug[k-1] + t * (n_aug[k] - n_aug[k-1])
            break
        end
    end
    isnan(δ99) && return (NaN, NaN, NaN)

    # δ*, θ trapezoid 0 → δ99
    δstar = 0.0; θmom = 0.0
    @inbounds for k in 2:length(n_aug)
        n1, n2 = n_aug[k-1], n_aug[k]
        r1 = u_aug[k-1] / Ue
        if n2 >= δ99
            t = (δ99 - n1) / (n2 - n1)
            r2 = r1 + t * (u_aug[k] / Ue - r1)
            dn = δ99 - n1
            δstar += 0.5 * ((1 - r1) + (1 - r2))       * dn
            θmom  += 0.5 * (r1*(1 - r1) + r2*(1 - r2)) * dn
            break
        else
            r2 = u_aug[k] / Ue
            dn = n2 - n1
            δstar += 0.5 * ((1 - r1) + (1 - r2))       * dn
            θmom  += 0.5 * (r1*(1 - r1) + r2*(1 - r2)) * dn
        end
    end
    return (δ99, δstar, θmom)
end

# ---------------------------------------------------------------------------
# Main (Pass 4: + per-face integrals → blMetrics.csv)
# ---------------------------------------------------------------------------

function main(case_dir::AbstractString)
    inputParam = joinpath(case_dir, "constant", "inputParam")
    midPlane   = joinpath(case_dir, "postProcessing", "midPlane.csv")
    wallCoords = joinpath(case_dir, "postProcessing", "wallCoordinates.csv")

    @printf "blMetrics: case = %s\n" abspath(case_dir)
    isfile(inputParam) || error("missing $(inputParam)")
    isfile(midPlane)   || error("missing $(midPlane) — run runPostProcess first")
    isfile(wallCoords) || error("missing $(wallCoords) — partial mode not active?")

    up = read_subdict(inputParam, "Upinlet")
    bl = read_subdict(inputParam, "blMetrics")
    method     = bl["method"]
    expHeightM = up["exportHeight"] * 1e-3                            # mm → m

    hdr, data = read_midplane(midPlane)
    col_omz   = findfirst(==("omz"), hdr)
    col_omz === nothing && error(
        "midPlane.csv has no `omz` column — re-run runPostProcess with the OF-side changes")
    xw, yw    = read_wall_coords(wallCoords)
    nW        = length(xw)

    # Wall pressure (kinematic) — same row order as wallCoordinates.csv since
    # both are populated by the same filtered loop in writeMidPlane.
    wallQ    = joinpath(case_dir, "postProcessing", "wallQuantities.csv")
    isfile(wallQ) || error("missing $(wallQ) — needed for pressureBernoulli")
    wq_data, _ = readdlm(wallQ, ',', Float64, '\n'; header=true)
    p_wall = wq_data[:, 3]                    # col 3 = p in (x,s,p,tau,yPlus)
    Uinf   = up["freeStreamVelocityStreamwise"]
    Uinf2  = Uinf * Uinf

    # Geometry + wall-face parameterisation
    geom    = build_geom(up)
    xis, ss = project_wall_to_xi(geom, xw, yw)
    perm    = sortperm(xis)                       # ascending xi order for output sort
    frames  = [surface_frame(geom, xis[j]) for j in 1:nW]

    # Voronoi: assign each mid-plane cell to its nearest wall face
    near    = assign_cells_to_walls(data[:,1], data[:,2], xw, yw)
    buckets = [Int[] for _ in 1:nW]
    @inbounds for i in eachindex(near); push!(buckets[near[i]], i); end

    # Per-face BL integrals — compute ALL three methods every run
    method in METHODS ||
        error("blMetrics.method = $(method) — supported: $(METHODS)")
    xs = @view data[:,1]; ys = @view data[:,2]
    us = @view data[:,4]; vs = @view data[:,5]
    ωs = @view data[:, col_omz]

    # Per-method arrays
    Ue_vT  = fill(NaN, nW); d99_vT = fill(NaN, nW); dst_vT = fill(NaN, nW); th_vT  = fill(NaN, nW)
    Ue_vM  = fill(NaN, nW); d99_vM = fill(NaN, nW); dst_vM = fill(NaN, nW); th_vM  = fill(NaN, nW)
    Ue_m   = fill(NaN, nW); d99_m  = fill(NaN, nW); dst_m  = fill(NaN, nW); th_m   = fill(NaN, nW)
    Ue_f   = fill(NaN, nW); d99_f  = fill(NaN, nW); dst_f  = fill(NaN, nW); th_f   = fill(NaN, nW)
    Ue_b   = fill(NaN, nW); d99_b  = fill(NaN, nW); dst_b  = fill(NaN, nW); th_b   = fill(NaN, nW)

    for j in 1:nW
        nxj, nyj, txj, tyj = frames[j]
        Pjx, Pjy = xw[j], yw[j]
        idx = buckets[j]; K = length(idx)
        n_col = Vector{Float64}(undef, K)
        u_col = Vector{Float64}(undef, K)
        ω_col = Vector{Float64}(undef, K)
        @inbounds for (k, i) in enumerate(idx)
            dx = xs[i] - Pjx; dy = ys[i] - Pjy
            n_col[k] = dx*nxj + dy*nyj
            u_col[k] = us[i]*txj + vs[i]*tyj
            ω_col[k] = ωs[i]
        end
        # Drop the wall-extrap row (n ≈ 0): its ω_z is the extrapolated value
        # and can be unphysical. Then sort ascending in n once (compute_ue and
        # compute_bl_integrals require this for the fixedHeight variant).
        keep   = findall(>(1e-9), n_col)
        nperm  = sortperm(view(n_col, keep))
        n_int  = n_col[keep][nperm]
        u_int  = u_col[keep][nperm]
        omz_int = ω_col[keep][nperm]

        Ue_vT[j]  = compute_ue(n_int, u_int, omz_int, "vorticityIntegralTrapezoidal")
        Ue_vM[j] = compute_ue(n_int, u_int, omz_int, "vorticityIntegralMidpoint")
        Ue_m[j]  = compute_ue(n_int, u_int, omz_int, "maxProfile")
        Ue_f[j]  = compute_ue(n_int, u_int, omz_int, "fixedHeight")
        # Bernoulli (Prandtl): p_e ≈ p_wall; recover U_e from total head.
        Ue_b[j]  = sqrt(max(0.0, Uinf2 - 2.0 * p_wall[j]))
        d99_vT[j],  dst_vT[j],  th_vT[j]  = compute_bl_integrals(n_int, u_int, Ue_vT[j];  max_n=expHeightM)
        d99_vM[j], dst_vM[j], th_vM[j] = compute_bl_integrals(n_int, u_int, Ue_vM[j]; max_n=expHeightM)
        d99_m[j],  dst_m[j],  th_m[j]  = compute_bl_integrals(n_int, u_int, Ue_m[j];  max_n=expHeightM)
        d99_f[j],  dst_f[j],  th_f[j]  = compute_bl_integrals(n_int, u_int, Ue_f[j];  max_n=expHeightM)
        d99_b[j],  dst_b[j],  th_b[j]  = compute_bl_integrals(n_int, u_int, Ue_b[j];  max_n=expHeightM)
    end

    # Production output (configured method) — sorted by xi
    sel = method == "vorticityIntegralTrapezoidal" ? (Ue_vT, d99_vT, dst_vT, th_vT) :
          method == "vorticityIntegralMidpoint"    ? (Ue_vM, d99_vM, dst_vM, th_vM) :
          method == "maxProfile"                   ? (Ue_m,  d99_m,  dst_m,  th_m) :
          method == "fixedHeight"                  ? (Ue_f,  d99_f,  dst_f,  th_f) :
                                                     (Ue_b,  d99_b,  dst_b,  th_b)
    Ue_sel, d99_sel, dst_sel, th_sel = sel

    outFile = joinpath(case_dir, "postProcessing", "BLQuantities.csv")
    open(outFile, "w") do io
        println(io, "x,s,xi,Ue,d99,dstar,Theta")
        for k in 1:nW
            j = perm[k]
            @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                    xw[j], ss[j], xis[j], Ue_sel[j], d99_sel[j], dst_sel[j], th_sel[j])
        end
    end

    # Comparison output (all five methods) — same row order
    cmpFile = joinpath(case_dir, "postProcessing", "BLQuantities_compare.csv")
    open(cmpFile, "w") do io
        println(io, "x,s," *
            "Ue_vorticityTrap,Ue_vorticityMid,Ue_max,Ue_fix,Ue_bern," *
            "d99_vorticityTrap,d99_vorticityMid,d99_max,d99_fix,d99_bern," *
            "dstar_vorticityTrap,dstar_vorticityMid,dstar_max,dstar_fix,dstar_bern," *
            "Theta_vorticityTrap,Theta_vorticityMid,Theta_max,Theta_fix,Theta_bern")
        for k in 1:nW
            j = perm[k]
            @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                xw[j], ss[j],
                Ue_vT[j],   Ue_vM[j],   Ue_m[j],   Ue_f[j],   Ue_b[j],
                d99_vT[j],  d99_vM[j],  d99_m[j],  d99_f[j],  d99_b[j],
                dst_vT[j],  dst_vM[j],  dst_m[j],  dst_f[j],  dst_b[j],
                th_vT[j],   th_vM[j],   th_m[j],   th_f[j],   th_b[j])
        end
    end

    # Summary — table across all five methods, the configured one is starred
    nW_valid(v) = count(isfinite, v)
    star(m) = (m == method ? " *" : "")
    @printf "  wall faces=%d  xi∈[%.4f, %.4f]   configured method = %s\n" nW minimum(xis) maximum(xis) method
    @printf "  %-26s   %-7s   %-8s   %-8s   %-8s   %-7s\n" "method" "valid" "median Ue" "med d99" "med d*" "med H"
    for (lbl, Ue_, d99_, ds_, th_) in (
            ("vorticityIntegralTrapezoidal", Ue_vT, d99_vT, dst_vT, th_vT),
            ("vorticityIntegralMidpoint",    Ue_vM, d99_vM, dst_vM, th_vM),
            ("maxProfile",                   Ue_m,  d99_m,  dst_m,  th_m),
            ("fixedHeight",                  Ue_f,  d99_f,  dst_f,  th_f),
            ("pressureBernoulli",            Ue_b,  d99_b,  dst_b,  th_b))
        ok = isfinite.(Ue_) .& isfinite.(d99_) .& isfinite.(ds_) .& isfinite.(th_) .& (Ue_ .> 0)
        nok = count(ok)
        if nok > 0
            @printf "  %-26s%s  %4d/%d  %8.3f   %7.4f   %7.4f   %5.3f\n" lbl star(lbl) nok nW median(Ue_[ok]) 1e3*median(d99_[ok]) 1e3*median(ds_[ok]) median(ds_[ok] ./ th_[ok])
        else
            @printf "  %-26s%s  %4d/%d  (no valid stations)\n" lbl star(lbl) 0 nW
        end
    end
    @printf "  → %s\n" outFile
    @printf "  → %s\n" cmpFile

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    case_dir = length(ARGS) >= 1 ? ARGS[1] : "."
    main(case_dir)
end
