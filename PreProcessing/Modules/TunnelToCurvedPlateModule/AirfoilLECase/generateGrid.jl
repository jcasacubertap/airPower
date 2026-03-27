#!/usr/bin/env julia
#=
    generate_cgrid.jl

    C-grid geometry generator for OpenFOAM blockMeshDict.
    Case: 221airfoilTopBottom (M3J symmetric airfoil, both surfaces)

    The C-grid wraps from the suction (upper) outlet, around the
    leading edge, to the pressure (lower) outlet.  Two radial layers
    (boundary-layer + far-field) give wall-normal resolution.

    Block topology (side view, one z-plane):

        outer ─────────────────────────────────────── outer
          │  far N  │ … │ far 2 │ far 1 │ far 0 │ … │ far M  │
        mid ─────────────────────────────────────── mid
          │  BL  N  │ … │ BL  2 │ BL  1 │ BL  0 │ … │ BL  M  │
        airfoil ══════════╗           ╔══════════ airfoil
        pressure outlet    ╚═══ LE ═══╝    suction inlet

    The C-path visits stations:
        suction outlet → xi_arch (upper) → LE → xi_arch (lower) → pressure outlet

    Usage:   julia generate_cgrid.jl [output_file]
    Default: writes  blockMeshDict_cgrid  in the current directory.
=#

using Printf

# ═══════════════════════════════════════════════════════════════════════
#  INPUTS — load inputs.jl (single source of truth)
# ═══════════════════════════════════════════════════════════════════════

"""Walk up from `start` until a file named `name` is found."""
function find_ancestor_file(start::String, name::String)
    dir = start
    while true
        candidate = joinpath(dir, name)
        isfile(candidate) && return candidate
        parent = dirname(dir)
        parent == dir && error("$name not found above $start")
        dir = parent
    end
end

include(find_ancestor_file(@__DIR__, "inputs.jl"))

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 1 — AIRFOIL DATA  (ξ = x/c,  η = y/c,  LE → TE)
# ═══════════════════════════════════════════════════════════════════════
# Read from .dat file (TE upper → LE → TE lower), extract upper surface,
# sort LE → TE.  Lower side is the mirror: η_lower(ξ) = −η_upper(ξ).

"""
    load_airfoil_upper(filepath) → (xi, eta)

Read an airfoil .dat file (tab/space-separated x/c, y/c columns, full
contour TE→LE→TE) and return the upper-surface points sorted LE→TE.
"""
function load_airfoil_upper(filepath::String)
    xi_all = Float64[]
    eta_all = Float64[]
    for line in eachline(filepath)
        line = strip(line)
        isempty(line) && continue
        parts = split(line)
        length(parts) >= 2 || continue
        x = tryparse(Float64, parts[1])
        y = tryparse(Float64, parts[2])
        (x === nothing || y === nothing) && continue
        push!(xi_all, x)
        push!(eta_all, y)
    end

    # Extract upper surface (eta ≥ 0), deduplicate TE (xi=1, eta=0)
    mask = eta_all .>= 0.0
    xi_up  = xi_all[mask]
    eta_up = eta_all[mask]

    # Sort by xi ascending (LE → TE)
    perm = sortperm(xi_up)
    xi_up  = xi_up[perm]
    eta_up = eta_up[perm]

    # Remove duplicate points (e.g. two (0,0) or two (1,0))
    keep = [true; [!(xi_up[i] == xi_up[i-1] && eta_up[i] == eta_up[i-1]) for i in 2:length(xi_up)]]
    return xi_up[keep], eta_up[keep]
end

# Resolve .dat path from inputs
const _AIRFOIL_DAT = let
    dir = find_ancestor_file(@__DIR__, "inputs.jl") |> dirname
    joinpath(dir, "PreProcessing", "InputOutput", "AirfoilGeometryData", inp.TTCP.tunnel.airfoilFile)
end
const XI_DATA, ETA_DATA = load_airfoil_upper(_AIRFOIL_DAT)

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 2 — LOAD USER PARAMETERS INTO GLOBALS
# ═══════════════════════════════════════════════════════════════════════

function load_params()
    t = inp.TTCP.tunnel
    a = inp.TTCP.airfoilLE
    mm(x) = x * 1000.0   # inputs.jl stores [m], generateGrid works in [mm]

    # — Airfoil characteristics —
    global CHORD      = mm(t.chord)
    global ALPHA_DEG  = Float64(t.alphaDeg)
    global X_CENTER   = mm(t.xCenter)
    global Y_CENTER   = mm(t.yCenter)

    # — Domain size —
    global XI_ARCH             = Float64(a.xiArch)
    global XI_SUCTION_OUTLET   = Float64(a.xiSuctionOutlet)
    global XI_PRESSURE_OUTLET  = Float64(a.xiPressureOutlet)
    global H_BL       = mm(a.hBL)
    global H_FAR      = mm(a.hFar)
    global Z_WIDTH    = 2.0           # [mm] spanwise extent (quasi-2D, cyclic)
    global NS_SUCTION  = Int(a.nSegSuction)
    global NS_ARCH_UP  = Int(a.nSegArchUp)
    global NS_ARCH_LO  = Int(a.nSegArchLo)
    global NS_PRESSURE = Int(a.nSegPressure)
    global COSINE_ARCH = Bool(a.cosineArch)

    # — Grid resolution —
    global NX_TOTAL    = Int(a.NxTotal)
    global NY_BL       = Int(a.NyBL)
    global NZ          = 2              # spanwise cells (quasi-2D)
    global GRAD_BL     = Float64(a.gradBL)
    global GRAD_ARCH   = Float64(a.gradArch)
end

# Declare globals (assigned by load_params)
CHORD = 0.0;  ALPHA_DEG = 0.0;  X_CENTER = 0.0;  Y_CENTER = 0.0
XI_ARCH = 0.0;  XI_SUCTION_OUTLET = 0.0;  XI_PRESSURE_OUTLET = 0.0
H_BL = 0.0;  H_FAR = 0.0;  Z_WIDTH = 0.0
NS_SUCTION = 0;  NS_ARCH_UP = 0;  NS_ARCH_LO = 0;  NS_PRESSURE = 0
COSINE_ARCH = true
NX_TOTAL = 0;  NY_BL = 0;  NZ = 0
GRAD_BL = 0.0;  GRAD_ARCH = 1.0
# Computed at runtime:
NY_FAR = 0;  GRAD_FAR = 1.0
NX_SEG = Int[]   # per-segment streamwise cell counts
GX_SEG = Float64[]  # per-segment streamwise grading


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 3 — INTERPOLATION & GEOMETRY HELPERS
# ═══════════════════════════════════════════════════════════════════════

"""Piecewise-linear interpolation in a sorted table."""
function interp1(xd::Vector{Float64}, yd::Vector{Float64}, xq::Float64)::Float64
    n = length(xd)
    xq ≤ xd[1]  && return yd[1]
    xq ≥ xd[n]  && return yd[n]
    for i in 1:n-1
        if xd[i] ≤ xq ≤ xd[i+1]
            t = (xq - xd[i]) / (xd[i+1] - xd[i])
            return (1.0 - t) * yd[i] + t * yd[i+1]
        end
    end
    return yd[n]
end

eta_upper(xi::Float64) =  interp1(XI_DATA, ETA_DATA, xi)
eta_lower(xi::Float64) = -interp1(XI_DATA, ETA_DATA, xi)

"""Cubic Lagrange interpolation using 4 nearest points from a sorted table."""
function interp1_cubic(xd::Vector{Float64}, yd::Vector{Float64}, xq::Float64)::Float64
    n = length(xd)
    xq ≤ xd[1]  && return yd[1]
    xq ≥ xd[n]  && return yd[n]
    k = 1
    for i in 1:n-1
        if xd[i] ≤ xq ≤ xd[i+1]; k = i; break; end
    end
    i0 = clamp(k - 1, 1, n - 3)
    val = 0.0
    for i in 0:3
        basis = 1.0
        for j in 0:3
            j != i && (basis *= (xq - xd[i0+j]) / (xd[i0+i] - xd[i0+j]))
        end
        val += yd[i0+i] * basis
    end
    return val
end

upper_pt_cubic(xi::Float64) = to_global(xi,  interp1_cubic(XI_DATA, ETA_DATA, xi))
lower_pt_cubic(xi::Float64) = to_global(xi, -interp1_cubic(XI_DATA, ETA_DATA, xi))

"""
    to_global(xi, eta) → (x, y)  [mm]

Scale by chord, rotate by AoA about mid-chord, translate to (X_CENTER,Y_CENTER).
"""
function to_global(xi::Float64, eta::Float64)
    α  = deg2rad(ALPHA_DEG)
    ca, sa = cos(α), sin(α)
    Xp = xi  * CHORD - 0.5 * CHORD        # relative to mid-chord
    Yp = eta * CHORD
    Xr =  ca * Xp + sa * Yp               # rotate by -α (nose up for +AoA)
    Yr = -sa * Xp + ca * Yp
    return (X_CENTER + Xr, Y_CENTER + Yr)
end

upper_pt(xi::Float64) = to_global(xi, eta_upper(xi))
lower_pt(xi::Float64) = to_global(xi, eta_lower(xi))

"""
    surface_normal(xi, side) → (nx, ny)

Outward unit normal on the airfoil.
`side` ∈ {:upper, :lower}.  Upper normals point away from the airfoil
centre (generally upward); lower normals point generally downward.
"""
function surface_normal(xi::Float64, side::Symbol; δ::Float64=1e-5)
    pt_fn = side == :upper ? upper_pt : lower_pt
    x0, y0 = pt_fn(max(0.0, xi - δ))
    x1, y1 = pt_fn(min(1.0, xi + δ))
    tx, ty = x1 - x0, y1 - y0
    # +90° for upper (outward = up), −90° for lower (outward = down)
    nx, ny = side == :upper ? (-ty, tx) : (ty, -tx)
    mag = sqrt(nx^2 + ny^2)
    mag < 1e-15 && return side == :upper ? (0.0, 1.0) : (0.0, -1.0)
    return (nx / mag, ny / mag)
end

"""Cubic-interpolation version of surface_normal."""
function surface_normal_cubic(xi::Float64, side::Symbol; δ::Float64=1e-5)
    pt_fn = side == :upper ? upper_pt_cubic : lower_pt_cubic
    x0, y0 = pt_fn(max(0.0, xi - δ))
    x1, y1 = pt_fn(min(1.0, xi + δ))
    tx, ty = x1 - x0, y1 - y0
    nx, ny = side == :upper ? (-ty, tx) : (ty, -tx)
    mag = sqrt(nx^2 + ny^2)
    mag < 1e-15 && return side == :upper ? (0.0, 1.0) : (0.0, -1.0)
    return (nx / mag, ny / mag)
end

"""Cosine distribution: maps uniform s ∈ [0,1] → clustered s ∈ [0,1],
   with higher density near s=0 and s=1 (i.e. near the endpoints)."""
cosine_cluster(s::Float64) = 0.5 * (1.0 - cos(π * s))

"""
    equal_arc_xis(xi_from, xi_to, n_seg, side) → Vector{Float64}

Compute `n_seg + 1` ξ values between `xi_from` and `xi_to` (inclusive)
such that the arc length along the airfoil surface is equal for every segment.
`side` ∈ {:upper, :lower}.  The returned ξ values go from `xi_from` to `xi_to`.
"""
function equal_arc_xis(xi_from::Float64, xi_to::Float64, n_seg::Int,
                       side::Symbol)::Vector{Float64}
    pt_fn = side == :upper ? upper_pt : lower_pt

    # 1. Gather M3J data points in [min, max], plus endpoints
    xi_lo, xi_hi = minmax(xi_from, xi_to)
    xis = Float64[]
    for xi_d in XI_DATA
        if xi_lo ≤ xi_d ≤ xi_hi
            push!(xis, xi_d)
        end
    end
    # Ensure endpoints are included
    if isempty(xis) || abs(xis[1] - xi_lo) > 1e-12
        push!(xis, xi_lo)
    end
    if abs(xis[end] - xi_hi) > 1e-12
        push!(xis, xi_hi)
    end
    sort!(unique!(xis))

    # If C-path direction is decreasing ξ (upper surface), reverse
    going_down = xi_from > xi_to
    if going_down
        reverse!(xis)
    end

    # 2. Cumulative arc length along C-path order
    np = length(xis)
    cum = zeros(np)
    for k in 2:np
        p0 = pt_fn(xis[k-1])
        p1 = pt_fn(xis[k])
        cum[k] = cum[k-1] + sqrt((p1[1]-p0[1])^2 + (p1[2]-p0[2])^2)
    end
    total_arc = cum[end]

    # 3. Interpolate to find ξ at each equal-arc boundary
    result = Float64[xi_from]
    for seg in 1:n_seg-1
        target = seg * total_arc / n_seg
        for j in 2:np
            if cum[j] ≥ target - 1e-12
                t = (target - cum[j-1]) / (cum[j] - cum[j-1])
                xi_interp = xis[j-1] + t * (xis[j] - xis[j-1])
                push!(result, xi_interp)
                break
            end
        end
    end
    push!(result, xi_to)
    return result
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 4 — C-PATH STATIONS
# ═══════════════════════════════════════════════════════════════════════

struct Station
    xi   :: Float64       # ξ = x/c on the airfoil
    side :: Symbol        # :upper, :lower, or :le
    px   :: Float64       # surface point x  [mm]
    py   :: Float64       # surface point y  [mm]
    nx   :: Float64       # outward unit normal x
    ny   :: Float64       # outward unit normal y
end

"""
Build the ordered list of C-path stations:
    suction outlet → xi_arch (upper) → LE → xi_arch (lower) → pressure outlet

Boundary stations are shared between adjacent sections so the LE
appears exactly once and no zero-length segments are created.

    Section A: NS_SUCTION + 1   stations  (suction outlet → xi_arch,  upper)
    Section B: NS_ARCH_UP − 1   stations  (interior between xi_arch and LE)
    Section C: 1 station                   (LE, averaged normal)
    Section D: NS_ARCH_LO − 1   stations  (interior between LE and xi_arch)
    Section E: NS_PRESSURE + 1  stations  (xi_arch → pressure outlet,  lower)

Total stations  Nv = NS_SUCTION + NS_ARCH_UP + NS_ARCH_LO + NS_PRESSURE + 1
Total segments  Ns = Nv − 1

Sections A and E use **equal arc-length** spacing so that all segments
within each section have identical arc lengths (→ zero cell-size mismatch
within sections).  Sections B and D keep cosine clustering for LE resolution.
"""
function build_stations()::Vector{Station}
    sts = Station[]

    # ---- A. Upper (suction) surface: XI_SUCTION_OUTLET → XI_ARCH ----
    #       Equal arc-length spacing → identical arc per segment
    suction_xis = equal_arc_xis(XI_SUCTION_OUTLET, XI_ARCH, NS_SUCTION, :upper)
    for xi in suction_xis
        px, py = upper_pt(xi)
        nx, ny = surface_normal(xi, :upper)
        push!(sts, Station(xi, :upper, px, py, nx, ny))
    end

    # ---- B. Upper arch INTERIOR: between XI_ARCH and LE  ----
    for i in 1:NS_ARCH_UP - 1                                       # skip endpoints
        s  = i / NS_ARCH_UP
        sc = COSINE_ARCH ? cosine_cluster(s) : s                    # optional clustering
        xi = XI_ARCH * (1.0 - sc)                                   # → 0  (exclusive)
        px, py = upper_pt(xi)
        nx, ny = surface_normal(max(xi, 1e-6), :upper)
        push!(sts, Station(xi, :upper, px, py, nx, ny))
    end

    # ---- C. Leading-edge station (averaged upper/lower normal → upstream) ----
    px, py = upper_pt(0.0)
    nu = surface_normal(1e-6, :upper)
    nl = surface_normal(1e-6, :lower)
    nx = 0.5 * (nu[1] + nl[1])
    ny = 0.5 * (nu[2] + nl[2])
    mag = sqrt(nx^2 + ny^2)
    if mag > 1e-15; nx /= mag; ny /= mag; end
    push!(sts, Station(0.0, :le, px, py, nx, ny))

    # ---- D. Lower arch INTERIOR: between LE and XI_ARCH  ----
    for i in 1:NS_ARCH_LO - 1                                       # skip endpoints
        s  = i / NS_ARCH_LO
        sc = COSINE_ARCH ? cosine_cluster(s) : s
        xi = sc * XI_ARCH                                            # 0 → XI_ARCH (exclusive)
        px, py = lower_pt(xi)
        nx, ny = surface_normal(max(xi, 1e-6), :lower)
        push!(sts, Station(xi, :lower, px, py, nx, ny))
    end

    # ---- E. Lower (pressure) surface: XI_ARCH → XI_PRESSURE_OUTLET ----
    #       Equal arc-length spacing → identical arc per segment
    pressure_xis = equal_arc_xis(XI_ARCH, XI_PRESSURE_OUTLET, NS_PRESSURE, :lower)
    for xi in pressure_xis
        px, py = lower_pt(xi)
        nx, ny = surface_normal(max(xi, 1e-6), :lower)
        push!(sts, Station(xi, :lower, px, py, nx, ny))
    end

    return sts
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 4b — CELL-SIZE MATCHING
# ═══════════════════════════════════════════════════════════════════════

"""
    segment_arc_length(sts, si, sj)

Compute the surface arc length between stations `si` and `sj` (1-based)
by summing Euclidean distances through all M3J data points in between.
"""
function segment_arc_length(sts::Vector{Station}, si::Int, sj::Int)::Float64
    a = sts[si]
    b = sts[sj]
    on_lower = (a.side == :lower || (a.side == :le && b.side == :lower))

    xi_lo = min(a.xi, b.xi)
    xi_hi = max(a.xi, b.xi)

    # Build ordered chain:  [station_a, M3J interior points..., station_b]
    chain = Tuple{Float64,Float64}[(a.px, a.py)]

    xi_interior = Float64[]
    for xi_d in XI_DATA
        if xi_lo < xi_d < xi_hi
            push!(xi_interior, xi_d)
        end
    end
    if !on_lower
        sort!(xi_interior, rev=true)
    else
        sort!(xi_interior)
    end

    for xi in xi_interior
        if on_lower
            push!(chain, lower_pt(xi))
        else
            push!(chain, upper_pt(xi))
        end
    end
    push!(chain, (b.px, b.py))

    arc = 0.0
    for k in 2:length(chain)
        dx = chain[k][1] - chain[k-1][1]
        dy = chain[k][2] - chain[k-1][2]
        arc += sqrt(dx^2 + dy^2)
    end
    return arc
end

"""
    solve_r_for_d0(arc, d0, n) → r

Find the common ratio `r` of a geometric series of `n` terms, with first
term `d0`, that sums to `arc`:  d0 × (rⁿ − 1)/(r − 1) = arc.
Uses Newton's method.  Returns 1.0 when n ≤ 1 or the series is uniform.
"""
function solve_r_for_d0(arc::Float64, d0::Float64, n::Int)::Float64
    n ≤ 1 && return 1.0
    β = arc / d0                             # target: (rⁿ-1)/(r-1) = β
    abs(β - n) / max(n, 1) < 1e-10 && return 1.0   # uniform case

    # Initial guess from linear expansion: S ≈ n + n(n-1)/2·(r-1)
    r = 1.0 + 2.0 * (β / n - 1.0) / max(n - 1, 1)
    r = clamp(r, 0.01, 100.0)

    for _ in 1:100
        if abs(r - 1.0) < 1e-12
            f  = n - β
            fp = n * (n - 1) / 2.0
        else
            rn = r^n
            s  = (rn - 1.0) / (r - 1.0)
            f  = s - β
            fp = (n * r^(n-1) * (r - 1.0) - (rn - 1.0)) / (r - 1.0)^2
        end
        abs(f) < 1e-12 && break
        abs(fp) < 1e-30 && break
        r -= f / fp
        r = clamp(r, 1e-6, 1000.0)
    end
    return r
end

"""
    compute_cell_distribution(sts) → (nx_seg, gx_seg)

Compute per-segment cell counts and streamwise grading ratios.

**Straight sections** (suction, pressure): uniform grading gx = 1.0.
Equal-arc station spacing guarantees identical cell sizes within each section.

**Arch sections** — two kinds of block:
  • *Transition blocks* (outermost arch block on each side, adjacent to the
    straight section):  **auto-computed** grading via `solve_r_for_d0`.
    The first cell exactly equals the straight-section cell (chain input)
    and the grading smoothly bridges to the graded arch blocks.
  • *Graded blocks* (inner arch blocks touching the LE):  graded with the
    full `gradArch` compression/expansion to cluster cells toward the LE.

When `nSegArch == 1` there is no room for a separate transition block, so the
single block receives the full grading.

The cell-size chain propagates from suction through the arch to pressure,
guaranteeing **exact interface matching** at every chained boundary.  The only
rounding mismatch occurs at the pressure-section entry where the uniform cell
count rounds to the nearest integer.
"""
function compute_cell_distribution(sts::Vector{Station})::Tuple{Vector{Int}, Vector{Float64}}
    Ns = length(sts) - 1
    arcs = [segment_arc_length(sts, i, i+1) for i in 1:Ns]

    # Segment index ranges (1-based)
    i_suc = 1 : NS_SUCTION
    i_aup = NS_SUCTION+1 : NS_SUCTION+NS_ARCH_UP
    i_alo = NS_SUCTION+NS_ARCH_UP+1 : NS_SUCTION+NS_ARCH_UP+NS_ARCH_LO
    i_prs = NS_SUCTION+NS_ARCH_UP+NS_ARCH_LO+1 : Ns

    arc_suc = arcs[i_suc[1]]           # all equal (equal-arc spacing)
    arc_prs = arcs[i_prs[1]]           # all equal

    # Which arch blocks are transition (uniform) vs graded?
    has_trans_up = NS_ARCH_UP >= 2
    has_trans_lo = NS_ARCH_LO >= 2
    i_aup_vec = collect(i_aup)
    i_alo_vec = collect(i_alo)
    i_grad_up = has_trans_up ? i_aup_vec[2:end] : i_aup_vec
    i_grad_lo = has_trans_lo ? i_alo_vec[1:end-1] : i_alo_vec
    total_grad_up_arc = sum(arcs[k] for k in i_grad_up)
    total_grad_lo_arc = sum(arcs[k] for k in i_grad_lo)

    dx_nom = sum(arcs) / NX_TOTAL

    best_score  = Inf
    best_nx_seg = zeros(Int, Ns)
    best_gx_seg = ones(Ns)

    # ---- 3-D search over (nx_suc, nt_up, nt_lo) ----
    search_half = max(30, round(Int, NX_TOTAL / 12))
    nx_suc_lo = max(1, round(Int, arc_suc / dx_nom) - search_half)
    nx_suc_hi = round(Int, arc_suc / dx_nom) + search_half

    for nx_suc in nx_suc_lo:nx_suc_hi
        dx_suc = arc_suc / nx_suc

        # Upper transition search range
        if has_trans_up
            arc_tu = arcs[i_aup_vec[1]]
            nt_up_nom = max(2, round(Int, arc_tu / dx_suc))
            nt_up_range = max(2, nt_up_nom - 5) : nt_up_nom + 5
        else
            nt_up_range = 0:0       # single dummy iteration
        end

        for nt_up in nt_up_range

            # ---- Build upper arch chain ----
            d0 = dx_suc
            nx_aup = Int[];  gx_aup = Float64[]

            for (j, k) in enumerate(i_aup_vec)
                arc_k = arcs[k]
                if j == 1 && has_trans_up
                    r  = solve_r_for_d0(arc_k, d0, nt_up)
                    gx = r^(nt_up - 1)
                    push!(nx_aup, nt_up);  push!(gx_aup, gx)
                    d0 *= gx
                else
                    frac   = clamp(arc_k / total_grad_up_arc, 0.0, 1.0)
                    tgt_gx = (1.0 / GRAD_ARCH)^frac
                    avg  = d0 * (1.0 + tgt_gx) / 2.0
                    nest = max(2, round(Int, arc_k / avg))
                    best_n = nest;  best_err = Inf
                    for n in max(2, nest-20) : nest+20
                        r  = solve_r_for_d0(arc_k, d0, n)
                        gx = r^(n - 1)
                        err = abs(log(max(gx, 1e-30)) - log(max(tgt_gx, 1e-30)))
                        if err < best_err;  best_err = err;  best_n = n;  end
                    end
                    r  = solve_r_for_d0(arc_k, d0, best_n)
                    gx = r^(best_n - 1)
                    push!(nx_aup, best_n);  push!(gx_aup, gx)
                    d0 *= gx
                end
            end
            d_LE = d0

            # ---- Lower arch GRADED blocks ----
            d0 = d_LE
            nx_alo = Int[];  gx_alo = Float64[]

            graded_lo = has_trans_lo ? i_alo_vec[1:end-1] : i_alo_vec
            for k in graded_lo
                arc_k = arcs[k]
                frac   = clamp(arc_k / total_grad_lo_arc, 0.0, 1.0)
                tgt_gx = GRAD_ARCH^frac
                avg  = d0 * (1.0 + tgt_gx) / 2.0
                nest = max(2, round(Int, arc_k / avg))
                best_n = nest;  best_err = Inf
                for n in max(2, nest-20) : nest+20
                    r  = solve_r_for_d0(arc_k, d0, n)
                    gx = r^(n - 1)
                    err = abs(log(max(gx, 1e-30)) - log(max(tgt_gx, 1e-30)))
                    if err < best_err;  best_err = err;  best_n = n;  end
                end
                r  = solve_r_for_d0(arc_k, d0, best_n)
                gx = r^(best_n - 1)
                push!(nx_alo, best_n);  push!(gx_alo, gx)
                d0 *= gx
            end
            d_graded_end = d0

            # ---- Lower transition + pressure: search nt_lo ----
            if has_trans_lo
                arc_trans_lo = arcs[i_alo_vec[end]]
                nt_lo_nom = max(2, round(Int, arc_trans_lo / d_graded_end))

                for nt_lo in max(2, nt_lo_nom - 15) : nt_lo_nom + 15
                    r_t  = solve_r_for_d0(arc_trans_lo, d_graded_end, nt_lo)
                    gx_t = r_t^(nt_lo - 1)
                    (gx_t < 0.3 || gx_t > 3.0) && continue
                    d_last_t = d_graded_end * gx_t

                    nx_prs = max(1, round(Int, arc_prs / d_last_t))
                    dx_prs = arc_prs / nx_prs

                    max_mis = abs(dx_prs - d_last_t) / max(d_last_t, 1e-30)
                    total = NS_SUCTION * nx_suc + sum(nx_aup) +
                            sum(nx_alo) + nt_lo + NS_PRESSURE * nx_prs
                    total_deviation = abs(total - NX_TOTAL) / NX_TOTAL
                    score = max_mis + 0.1 * total_deviation

                    if score < best_score
                        best_score = score
                        fill!(best_nx_seg, 0);  fill!(best_gx_seg, 1.0)
                        for i in i_suc;  best_nx_seg[i] = nx_suc;  end
                        for (j, i) in enumerate(i_aup_vec)
                            best_nx_seg[i] = nx_aup[j];  best_gx_seg[i] = gx_aup[j]
                        end
                        for (j, kk) in enumerate(graded_lo)
                            best_nx_seg[kk] = nx_alo[j];  best_gx_seg[kk] = gx_alo[j]
                        end
                        best_nx_seg[i_alo_vec[end]] = nt_lo
                        best_gx_seg[i_alo_vec[end]] = gx_t
                        for i in i_prs;  best_nx_seg[i] = nx_prs;  end
                    end
                end
            else
                nx_prs = max(1, round(Int, arc_prs / d_graded_end))
                dx_prs = arc_prs / nx_prs

                max_mis = abs(dx_prs - d_graded_end) / max(d_graded_end, 1e-30)
                total = NS_SUCTION * nx_suc + sum(nx_aup) +
                        sum(nx_alo) + NS_PRESSURE * nx_prs
                total_deviation = abs(total - NX_TOTAL) / NX_TOTAL
                score = max_mis + 0.1 * total_deviation

                if score < best_score
                    best_score = score
                    fill!(best_nx_seg, 0);  fill!(best_gx_seg, 1.0)
                    for i in i_suc;  best_nx_seg[i] = nx_suc;  end
                    for (j, i) in enumerate(i_aup_vec)
                        best_nx_seg[i] = nx_aup[j];  best_gx_seg[i] = gx_aup[j]
                    end
                    for (j, i) in enumerate(i_alo_vec)
                        best_nx_seg[i] = nx_alo[j];  best_gx_seg[i] = gx_alo[j]
                    end
                    for i in i_prs;  best_nx_seg[i] = nx_prs;  end
                end
            end

        end  # nt_up
    end  # nx_suc

    return best_nx_seg, best_gx_seg
end

"""
    compute_ny_far()

Compute NY_FAR so the first far-field cell matches the last BL cell.
With geometric grading in the BL block (expansion ratio GRAD_BL over NY_BL cells):
    common ratio   r = GRAD_BL^(1/(NY_BL-1))
    first cell     d₀ = H_BL × (r-1) / (r^NY_BL - 1)
    last BL cell   d_outer = d₀ × GRAD_BL

For a uniform far-field: NY_FAR = round(H_FAR / d_outer)
"""
function compute_ny_far()
    if GRAD_BL ≈ 1.0
        d_outer = H_BL / NY_BL
    else
        r = GRAD_BL^(1.0 / (NY_BL - 1))
        d0 = H_BL * (r - 1.0) / (r^NY_BL - 1.0)
        d_outer = d0 * GRAD_BL
    end
    ny = max(1, round(Int, H_FAR / d_outer))
    return ny, d_outer
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 5 — VERTICES
# ═══════════════════════════════════════════════════════════════════════
# Vertex layout  (Nv = number of stations):
#
#   Front z-plane  (z = −Z_WIDTH/2):
#       level 0  (surface):  indices  0  …  Nv−1
#       level 1  (mid/BL) :  indices  Nv …  2Nv−1
#       level 2  (outer)  :  indices  2Nv…  3Nv−1
#
#   Back z-plane   (z = +Z_WIDTH/2):
#       level 0 :  3Nv …  4Nv−1
#       level 1 :  4Nv …  5Nv−1
#       level 2 :  5Nv …  6Nv−1

struct Vertex
    x :: Float64
    y :: Float64
    z :: Float64
end

function compute_vertices(sts::Vector{Station})::Vector{Vertex}
    Nv   = length(sts)
    verts = Vertex[]
    sizehint!(verts, 6Nv)

    for (ip, zval) in enumerate([-Z_WIDTH/2, Z_WIDTH/2])       # front=1, back=2
        for H in [0.0, H_BL, H_BL + H_FAR]                    # surface, mid, outer
            for st in sts
                push!(verts, Vertex(st.px + H * st.nx,
                                    st.py + H * st.ny,
                                    zval))
            end
        end
    end
    return verts
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 6 — BLOCKS
# ═══════════════════════════════════════════════════════════════════════
# Each C-path segment  i → i+1  produces two hex blocks:
#   • near-wall  (level 0 → level 1)
#   • far-field  (level 1 → level 2)

struct HexBlock
    v   :: NTuple{8,Int}       # 0-based vertex indices
    nx  :: Int                 # cells in ξ
    ny  :: Int                 # cells in η
    nz  :: Int                 # cells in z
    gx  :: Float64             # grading in ξ
    gy  :: Float64             # grading in η
    gz  :: Float64             # grading in z
end

function compute_blocks(Nv::Int, nx_seg::Vector{Int},
                        gx_seg::Vector{Float64})::Vector{HexBlock}
    Ns = Nv - 1
    blks = HexBlock[]
    sizehint!(blks, 2Ns)

    for i in 0:Ns-1
        nxi = nx_seg[i+1]
        gxi = gx_seg[i+1]

        # Near-wall block  (back z-plane first → right-hand normal toward front)
        push!(blks, HexBlock(
            (3Nv+i, 3Nv+i+1, 4Nv+i+1, 4Nv+i,
             i,     i+1,     Nv+i+1,  Nv+i),
            nxi, NY_BL, NZ, gxi, GRAD_BL, 1.0))

        # Far-field block
        push!(blks, HexBlock(
            (4Nv+i, 4Nv+i+1, 5Nv+i+1, 5Nv+i,
             Nv+i,  Nv+i+1,  2Nv+i+1, 2Nv+i),
            nxi, NY_FAR, NZ, gxi, GRAD_FAR, 1.0))
    end
    return blks
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 7 — SPLINE EDGES
# ═══════════════════════════════════════════════════════════════════════
# For each segment, three horizontal edges per z-plane (surface, mid, outer).

struct SplineEdge
    v0   :: Int                            # start vertex (0-based)
    v1   :: Int                            # end   vertex (0-based)
    pts  :: Vector{NTuple{3,Float64}}      # intermediate points
    etype :: String                        # "spline" or "polyLine"
end

"""
Return `n-1` interior ξ values equally spaced by arc length between
`xi_from` and `xi_to` on the given surface side.  When `H > 0` the arc
length is measured along the **offset curve** (surface + H·normal) so that
the resulting spline knots give uniform cells at the actual radial height,
not just at the wall.
"""
function resample_arc(xi_from::Float64, xi_to::Float64, n::Int,
                      side::Symbol, H::Float64=0.0)::Vector{Float64}
    pt_fn = side == :upper ? upper_pt : lower_pt
    xi_lo, xi_hi = minmax(xi_from, xi_to)

    # Gather M3J data points in [xi_lo, xi_hi] plus endpoints
    xis = Float64[]
    for xi_d in XI_DATA
        if xi_lo ≤ xi_d ≤ xi_hi
            push!(xis, xi_d)
        end
    end
    if isempty(xis) || abs(xis[1] - xi_lo) > 1e-12
        push!(xis, xi_lo)
    end
    if isempty(xis) || abs(xis[end] - xi_hi) > 1e-12
        push!(xis, xi_hi)
    end
    sort!(unique!(xis))

    # If C-path direction is decreasing ξ (upper surface), reverse
    going_down = xi_from > xi_to
    if going_down
        reverse!(xis)
    end

    # Compute physical points at offset height H
    np = length(xis)
    ox = zeros(np)
    oy = zeros(np)
    for k in 1:np
        px, py = pt_fn(xis[k])
        if H > 0.0
            nx, ny = surface_normal(max(xis[k], 1e-6), side)
            ox[k] = px + H * nx
            oy[k] = py + H * ny
        else
            ox[k] = px
            oy[k] = py
        end
    end

    # Cumulative arc length along the offset curve
    cum = zeros(np)
    for k in 2:np
        cum[k] = cum[k-1] + sqrt((ox[k]-ox[k-1])^2 + (oy[k]-oy[k-1])^2)
    end
    total_arc = cum[end]

    # Interpolate n-1 interior points at equal arc-length intervals
    result = Float64[]
    for i in 1:n-1
        target = i * total_arc / n
        for j in 2:np
            if cum[j] ≥ target - 1e-12
                t = (target - cum[j-1]) / (cum[j] - cum[j-1])
                push!(result, xis[j-1] + t * (xis[j] - xis[j-1]))
                break
            end
        end
    end
    return result
end

# Arch blocks use polyLine (piecewise-linear edges) which need many
# knots to trace the curved surface.  Straight blocks use Catmull-Rom
# splines where fewer knots work well (too many causes overshoot).
const N_RESAMPLE_ARCH = max(50, round(Int, NX_TOTAL / 8))
const N_RESAMPLE_STRAIGHT = 10

"""
Compute intermediate spline points between two adjacent stations
at a given radial offset `H` and z-coordinate `zval`.

Interior ξ values are resampled at equal arc-length intervals so that
OpenFOAM's Catmull-Rom parameterization produces uniform cell sizes
across block interfaces.  Arch blocks keep original M3J data points
to preserve LE curvature.
"""
function spline_between(sts::Vector{Station}, si::Int, sj::Int,
                        H::Float64, zval::Float64)::Vector{NTuple{3,Float64}}
    a = sts[si]   # Julia 1-based
    b = sts[sj]

    # Determine which surface side we are on
    on_lower = (a.side == :lower || (a.side == :le && b.side == :lower))

    # Arch blocks at the wall (H≈0) keep the original M3J data points to
    # preserve the accurate LE curvature.  At outer radial levels (H>0),
    # all blocks use arc-length resampling at the offset height — this
    # avoids the fragile dependency on which M3J points happen to fall
    # inside each block (which changes with XI_ARCH and other parameters).
    is_arch_wall = (max(a.xi, b.xi) <= XI_ARCH) && (H < 1e-12)

    xi_lo = min(a.xi, b.xi)
    xi_hi = max(a.xi, b.xi)

    if is_arch_wall
        # Original M3J data points only — cubic surface evaluation (below)
        # provides the LE curvature; extra midpoints would constrain the
        # Catmull-Rom to the piecewise-linear shape and reduce smoothness.
        xi_interior = Float64[]
        for xi_d in XI_DATA
            if xi_lo < xi_d < xi_hi
                push!(xi_interior, xi_d)
            end
        end
    else
        # Arc-length resampling at offset height H for uniform knot spacing
        # Arch blocks need dense knots (polyLine); straight blocks need few (spline)
        in_arch_h = max(a.xi, b.xi) <= XI_ARCH
        n_base = in_arch_h ? N_RESAMPLE_ARCH : N_RESAMPLE_STRAIGHT
        n_pts = in_arch_h ? 2 * n_base : n_base
        xi_interior = resample_arc(
            a.xi, b.xi,
            n_pts,
            on_lower ? :lower : :upper,
            H
        )
    end

    # Order the interior points along the C-path direction:
    #   upper surface: decreasing ξ  (a.xi > b.xi)
    #   lower surface: increasing ξ  (a.xi < b.xi)
    if !on_lower
        sort!(xi_interior, rev=true)
    else
        sort!(xi_interior)
    end

    # --- LE densification: add cubic-interpolated knots near ξ=0 ---
    # The first M3J segment (ξ=0 to ξ₁=0.001313) is large relative to the
    # nose curvature.  Adding denser knots via cubic Lagrange interpolation
    # through the first 4 real M3J data points steepens the Catmull-Rom
    # tangent at the LE vertex, shrinking the nose kink angle.
    le_touch = (a.side == :le) ? 1 : (b.side == :le) ? 2 : 0
    if le_touch > 0
        xd = XI_DATA[1:4]
        yd = ETA_DATA[1:4]
        # Add 3 knots in the first M3J segment via cubic Lagrange
        for frac in [0.03, 0.1, 0.3, 0.6]
            xi_sub = xd[1] + frac * (xd[2] - xd[1])
            # Skip if already covered by an existing interior point
            if xi_lo < xi_sub < xi_hi && !any(abs(xi_sub - x) < 1e-10 for x in xi_interior)
                eta_sub = 0.0
                for i in 1:4
                    basis = 1.0
                    for j in 1:4
                        j != i && (basis *= (xi_sub - xd[j]) / (xd[i] - xd[j]))
                    end
                    eta_sub += yd[i] * basis
                end
                push!(xi_interior, xi_sub)
            end
        end
        # Re-sort after adding densified knots
        if !on_lower
            sort!(xi_interior, rev=true)
        else
            sort!(xi_interior)
        end
    end

    # Use cubic surface evaluation only for LE-touching blocks where the
    # nose curvature is highest.  At xi=0 the cubic=linear (data point),
    # so the mismatch is only at the other vertex (~0.06mm at xi≈0.01).
    # Always use LINEAR normals — cubic normals oscillate too much
    # (up to ±5° / 24mm at H=250) due to Lagrange stencil shifts.
    in_arch = max(a.xi, b.xi) <= XI_ARCH

    pts = NTuple{3,Float64}[]
    for xi in xi_interior
        if on_lower
            nx, ny = surface_normal(max(xi, 1e-6), :lower)
            px, py = in_arch ? lower_pt_cubic(xi) : lower_pt(xi)
        else
            nx, ny = surface_normal(max(xi, 1e-6), :upper)
            px, py = in_arch ? upper_pt_cubic(xi) : upper_pt(xi)
        end
        push!(pts, (px + H * nx, py + H * ny, zval))
    end
    return pts
end

function compute_edges(sts::Vector{Station})::Vector{SplineEdge}
    Nv  = length(sts)
    Ns  = Nv - 1
    edges = SplineEdge[]

    for (ip, zval) in enumerate([-Z_WIDTH/2, Z_WIDTH/2])
        z_offset = (ip - 1) * 3Nv              # 0 for front, 3Nv for back

        for (lev, H) in enumerate([0.0, H_BL, H_BL + H_FAR])
            lev_offset = (lev - 1) * Nv         # 0, Nv, 2Nv

            for i in 0:Ns-1
                vi = z_offset + lev_offset + i
                vj = z_offset + lev_offset + i + 1
                pts = spline_between(sts, i+1, i+2, H, zval)   # 1-based station indices
                # Use polyLine for all blocks — traces knot positions exactly
                # without Catmull-Rom overshoot at high resolution, and avoids
                # tangent mismatches at polyLine/spline block boundaries.
                etype = "polyLine"
                push!(edges, SplineEdge(vi, vj, pts, etype))
            end
        end
    end
    return edges
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 8 — PATCHES  (boundary face definitions)
# ═══════════════════════════════════════════════════════════════════════
# Face = 4 vertex indices (0-based), ordered for outward normal.

const Face = NTuple{4,Int}

struct PatchDef
    name  :: String
    ptype :: String        # "wall", "patch", etc.
    faces :: Vector{Face}
end

function compute_patches(Nv::Int)::Vector{PatchDef}
    Ns = Nv - 1
    patches = PatchDef[]

    # ---- airfoil (wall): inner C-surface, all segments, level 0 ----
    af = Face[]
    for i in 0:Ns-1
        push!(af, (i, 3Nv+i, 3Nv+i+1, i+1))
    end
    push!(patches, PatchDef("airfoil", "wall", af))

    # ---- farfield (outer C-boundary): level 2, all segments ----
    ff = Face[]
    for i in 0:Ns-1
        push!(ff, (2Nv+i, 2Nv+i+1, 5Nv+i+1, 5Nv+i))
    end
    push!(patches, PatchDef("farfield", "patch", ff))

    # ---- suctionOutlet: upper end of C (station 0), levels 0→2 ----
    si = Face[]
    push!(si, (0, Nv, 4Nv, 3Nv))
    push!(si, (Nv, 2Nv, 5Nv, 4Nv))
    push!(patches, PatchDef("suctionOutlet", "patch", si))

    # ---- pressureOutlet: lower end of C (station Nv-1), levels 0→2 ----
    po = Face[]
    j = Nv - 1
    push!(po, (j, j+Nv, j+4Nv, j+3Nv))
    push!(po, (j+Nv, j+2Nv, j+5Nv, j+4Nv))
    push!(patches, PatchDef("pressureOutlet", "patch", po))

    # ---- front: all block faces at z = −Z_WIDTH/2 ----
    fr = Face[]
    for i in 0:Ns-1
        push!(fr, (i, i+1, Nv+i+1, Nv+i))
        push!(fr, (Nv+i, Nv+i+1, 2Nv+i+1, 2Nv+i))
    end
    push!(patches, PatchDef("front", "empty", fr))

    # ---- back: all block faces at z = +Z_WIDTH/2 ----
    bk = Face[]
    for i in 0:Ns-1
        push!(bk, (3Nv+i, 3Nv+i+1, 4Nv+i+1, 4Nv+i))
        push!(bk, (4Nv+i, 4Nv+i+1, 5Nv+i+1, 5Nv+i))
    end
    push!(patches, PatchDef("back", "empty", bk))

    return patches
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 9 — WRITE blockMeshDict
# ═══════════════════════════════════════════════════════════════════════

function write_blockmeshdict(fname::String,
                             verts::Vector{Vertex},
                             blks::Vector{HexBlock},
                             edges::Vector{SplineEdge},
                             patches::Vector{PatchDef})
    open(fname, "w") do io
        # ---- Header ----
        println(io, """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Generated by generate_cgrid.jl  (221airfoilTopBottom)
//
// C-grid parameters:
//   chord      = $CHORD mm,  AoA = $ALPHA_DEG deg
//   xi_arch    = $XI_ARCH
//   xi_suction = $XI_SUCTION_OUTLET,  xi_pressure = $XI_PRESSURE_OUTLET
//   H_BL       = $H_BL mm,  H_far = $H_FAR mm
//   segments   = $(NS_SUCTION)+$(NS_ARCH_UP)+1+$(NS_ARCH_LO)+$(NS_PRESSURE)  (suction+archUp+LE+archLo+pressure)
//   total blocks = $(length(blks))
//   cosine arch = $COSINE_ARCH

scale 0.001;   // vertices in mm → metres
""")

        # ---- Vertices ----
        println(io, "vertices\n(")
        for (k, v) in enumerate(verts)
            @printf(io, "    ( %14.8f  %14.8f  %14.8f )   // %d\n",
                    v.x, v.y, v.z, k - 1)
        end
        println(io, ");\n")

        # ---- Blocks ----
        println(io, "blocks\n(")
        for (k, b) in enumerate(blks)
            vs = join(string.(b.v), " ")
            @printf(io, "    hex (%s)  (%d %d %d) simpleGrading (%g %g %g)   // block %d\n",
                    vs, b.nx, b.ny, b.nz, b.gx, b.gy, b.gz, k - 1)
        end
        println(io, ");\n")

        # ---- Edges ----
        println(io, "edges\n(")
        for e in edges
            print(io, "    $(e.etype) $(e.v0) $(e.v1) (")
            for p in e.pts
                @printf(io, " ( %14.8f %14.8f %14.8f )", p[1], p[2], p[3])
            end
            println(io, " )")
        end
        println(io, ");\n")

        # ---- Boundary ----
        println(io, "boundary\n(")
        for p in patches
            println(io, "    $(p.name)\n    {\n        type $(p.ptype);")
            println(io, "        faces\n        (")
            for f in p.faces
                println(io, "            ($(f[1]) $(f[2]) $(f[3]) $(f[4]))")
            end
            println(io, "        );\n    }\n")
        end
        println(io, ");\n")

        println(io, "mergePatchPairs\n(\n);\n")
        println(io, "// ********************************************************** //")
    end
end


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 10 — SUMMARY  (print station table for verification)
# ═══════════════════════════════════════════════════════════════════════

function print_summary(sts::Vector{Station}, Nv::Int, nx_seg::Vector{Int},
                       gx_seg::Vector{Float64}, ny_far::Int, d_outer::Float64)
    Ns = Nv - 1
    arcs = [segment_arc_length(sts, i, i+1) for i in 1:Ns]

    # Compute first/last cell sizes per segment (accounting for grading)
    d_first = zeros(Ns)
    d_last  = zeros(Ns)
    for i in 1:Ns
        n  = nx_seg[i]
        gx = gx_seg[i]
        if n ≤ 1 || abs(gx - 1.0) < 1e-12
            d_first[i] = arcs[i] / max(n, 1)
            d_last[i]  = d_first[i]
        else
            r = gx^(1.0 / (n - 1))
            d_first[i] = arcs[i] * (r - 1.0) / (r^n - 1.0)
            d_last[i]  = d_first[i] * gx
        end
    end

    println("\n===  C-grid summary  (221airfoilTopBottom)  ===")
    println("  Stations (Nv)       = $Nv")
    println("  Segments (Ns)       = $Ns")
    println("  Blocks              = $(2Ns)   ($(Ns) BL + $(Ns) far-field)")
    println("  Vertices            = $(6Nv)")
    println("  Spline edges        = $(6Ns)")
    println("  NX total            = $(sum(nx_seg))   (requested $NX_TOTAL)")
    println("  NY_BL               = $NY_BL   (gradBL = $GRAD_BL)")
    println("  NY_FAR              = $ny_far   (uniform, auto-matched)")
    println("  gradArch            = $GRAD_ARCH")
    @printf("  BL outer cell size  = %.4f mm\n", d_outer)
    @printf("  Far-field cell size = %.4f mm\n", H_FAR / ny_far)

    # Report interface mismatches (last cell of seg i vs first cell of seg i+1)
    if Ns > 1
        mis_abs = [abs(d_last[i] - d_first[i+1]) for i in 1:Ns-1]
        mis_rel = [mis_abs[i] / min(d_last[i], d_first[i+1]) for i in 1:Ns-1]
        worst = argmax(mis_abs)
        @printf("  Max interface Δ(dx) = %.6f mm  (%.4f%%)  at seg %d→%d\n",
                mis_abs[worst], 100*mis_rel[worst], worst-1, worst)
    end
    println()

    @printf("  %4s  %6s  %8s  %10s  %5s  %8s  %10s  %10s\n",
            "seg", "side", "ξ_start", "arc [mm]", "Nx", "gx", "d_first", "d_last")
    println("  ", "-"^80)
    for i in 1:Ns
        st = sts[i]
        @printf("  %4d  %6s  %8.5f  %10.4f  %5d  %8.4f  %10.4f  %10.4f\n",
                i-1, st.side, st.xi, arcs[i], nx_seg[i], gx_seg[i],
                d_first[i], d_last[i])
    end
    println()
end


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

function main()
    # 0. Read parameters from system/inputDomain
    load_params()

    outfile = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "system", "blockMeshDict")

    # 1. Build C-path stations
    sts = build_stations()
    Nv  = length(sts)

    # 2. Cell-size matching
    #    a) Streamwise: uniform in straight sections, graded in arch (LE clustering)
    nx_seg, gx_seg = compute_cell_distribution(sts)
    global NX_SEG = nx_seg
    global GX_SEG = gx_seg

    #    b) Wall-normal: match BL outer cell size → uniform far-field
    ny_far, d_outer = compute_ny_far()
    global NY_FAR  = ny_far
    global GRAD_FAR = 1.0

    # 3. Compute mesh entities
    verts   = compute_vertices(sts)
    blks    = compute_blocks(Nv, NX_SEG, GX_SEG)
    edges   = compute_edges(sts)
    patches = compute_patches(Nv)

    # 4. Print summary
    print_summary(sts, Nv, NX_SEG, GX_SEG, ny_far, d_outer)

    # 5. Write blockMeshDict
    write_blockmeshdict(outfile, verts, blks, edges, patches)
    println("Wrote $outfile  ($(length(verts)) vertices, $(length(blks)) blocks, $(length(edges)) splines)")
end

main()
