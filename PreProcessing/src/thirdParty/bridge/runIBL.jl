using MAT

"""
    _find_matlab() → String | nothing

Locate a MATLAB executable. Tries `PATH` first, then the standard macOS
application bundles (newest version wins).
"""
function _find_matlab()
    exe = Sys.which("matlab")
    exe !== nothing && return exe
    candidates = sort(filter(isfile,
        [joinpath(d, "bin", "matlab")
         for d in glob_app_dirs()]); rev=true)
    return isempty(candidates) ? nothing : candidates[1]
end

glob_app_dirs() = isdir("/Applications") ?
    [joinpath("/Applications", d) for d in readdir("/Applications")
     if startswith(d, "MATLAB") && endswith(d, ".app")] : String[]

"""
    _ue_poly(x, p) → Float64

Edge (freestream) velocity at physical arc-length `x` from the virtual leading
edge, evaluated from the SAME log-pressure polynomial the DFP imposes on its
top boundary and coded inlet BC (`0/U`, Casacuberta et al. 2022):

    Ue(x) = Uinf · P(L),   L = log(x),
    P(L)  = (((pa4·L + pa3)·L + pa2)·L + pa1)·L + pa0.

Using this (rather than the simulated top-row u) guarantees the IBL and the DFP
share the same Ue(x), so the Falkner-Skan-Cooke inlet parameter beta — and
hence the imposed inflow profile — coincide at x = 0.
"""
function _ue_poly(x::Float64, p)
    L = log(x)
    return p.Uinf * ((((p.pa4*L + p.pa3)*L + p.pa2)*L + p.pa1)*L + p.pa0)
end

"""
    run_ibl_solver() → NamedTuple | nothing

Run the spectral integral boundary-layer solver (solver_IBL_spectral, via
MATLAB) for the DirectFlatPlate case and return its velocity fields (in memory)
for overlay on the DFP plots. Returns `nothing` if the solver cannot be run
(e.g. MATLAB unavailable) so that visualization can proceed without the overlay.
Leaves no files behind (see `run_ibl`).

The spectral solver seeds u, v AND w at the inlet from the Falkner-Skan-Cooke
similarity solution, so the crossflow is imposed at x = 0 and matches the DFP
at every station (the Thomas-algorithm solver_IBL grows w from a uniform inflow
and does not).

Solver inputs are derived from `inp.DFP` and `inp.VAL`:
  S  = xInlet,  L = xInlet + domainLength,  H = domainHeight,
  nu = freeStreamViscosity,  We = Winf,
  nx  = inp.VAL.blnx (streamwise marching points),
  ny  = inp.VAL.blNcheb (number of wall-normal Chebyshev nodes),
  y_i = inp.VAL.blYi (Chebyshev node median — half the nodes below it),
  Ue  = analytic edge velocity `Uinf·P(log x)` — the SAME distribution the DFP
        imposes (see `_ue_poly`), so the inlet FSC beta coincides with the DFP.

Returned NamedTuple fields:
  Xgrid  streamwise coordinate in the DFP frame (0 -> domainLength)
  Y      wall-normal Chebyshev coordinate [m] (H -> 0)
  u,v,w  ny x nxout IBL velocity components (row 1 = freestream, row ny = wall)
"""
function run_ibl_solver()
    p = inp.DFP

    # Shared IBL grid (inp.VAL): blnx streamwise points, blNcheb Chebyshev
    # wall-normal nodes over [0, H] with median at blYi. H is this case's domain
    # height (DFP: domainHeight).
    nx  = Int(inp.VAL.blnx)
    ny  = Int(inp.VAL.blNcheb)
    y_i = Float64(inp.VAL.blYi)

    # Edge velocity: the SAME analytic distribution the DFP imposes (log-
    # pressure polynomial), sampled on the uniform streamwise grid x_OF in
    # [0, domainLength] → physical arc-length x = xInlet .. xInlet+domainLength.
    # The virtual origin S = xInlet is given directly (no calibration).
    xg = range(0.0, p.domainLength, length=nx)
    Ue = [_ue_poly(p.xInlet + xq, p) for xq in xg]

    @info "DFP IBL grid & edge velocity" nx=nx ny_cheb=ny y_i=y_i H=p.domainHeight Ue_inlet=round(Ue[1], digits=4) Ue_outlet=round(Ue[end], digits=4)

    return run_ibl(; S=p.xInlet, Lspan=p.domainLength, H=p.domainHeight,
                   nu=p.freeStreamViscosity, We=p.Winf, Ue=Ue,
                   nx=nx, ny=ny, y_i=y_i)
end

"""
    run_ibl(; S, Lspan, H, nu, We, Ue, nx, ny, y_i, calibrate_dstar=0.0)
        → NamedTuple | nothing

Low-level spectral-IBL run shared by the DFP and TTCP paths. The MATLAB
input/output `.mat` files are written to a temporary directory and removed
afterwards — the run leaves no trace; only the returned in-memory solution
survives (or `nothing` on any failure). The virtual origin is `S` directly,
unless `calibrate_dstar > 0`, in which case the driver picks `S` so the FSC
inflow δ* equals that target. `Ue` is the length-`nx` edge velocity sampled over
the physical span `Lspan`.
"""
function run_ibl(; S::Real, Lspan::Real, H::Real, nu::Real, We::Real,
                 Ue::AbstractVector, nx::Integer, ny::Integer, y_i::Real,
                 calibrate_dstar::Real=0.0)
    matlab = _find_matlab()
    if matlab === nothing
        @warn "MATLAB executable not found — skipping IBL (set inp.VAL.valBL=false to silence)."
        return nothing
    end
    length(Ue) == nx || error("run_ibl: Ue has $(length(Ue)) entries but nx=$nx")

    driver_dir = @__DIR__                          # PreProcessing/src/thirdParty/bridge
    tmp = mktempdir()
    try
        infile  = joinpath(tmp, "ibl_in.mat")
        outfile = joinpath(tmp, "ibl_out.mat")
        matwrite(infile, Dict(
            "nx"  => Float64(nx),     "ny"    => Float64(ny),
            "S"   => Float64(S),      "Lspan" => Float64(Lspan), "H" => Float64(H),
            "y_i" => Float64(y_i),
            "nu"  => Float64(nu),     "We"    => Float64(We),
            "Ue"  => Vector{Float64}(Ue),
            "calibrate_dstar" => Float64(calibrate_dstar),
        ))

        cmd = Cmd(`$matlab -batch "ibl_driver('$infile','$outfile')"`; dir=driver_dir)
        @info "Invoking MATLAB IBL driver..." matlab=matlab calibrate_dstar=calibrate_dstar
        try
            run(cmd)
        catch e
            @warn "MATLAB IBL driver failed — skipping." exception=e
            return nothing
        end
        if !isfile(outfile)
            @warn "IBL driver produced no output — skipping."
            return nothing
        end
        bl = load_ibl_solution(outfile)
        @info "IBL solution loaded" stations=length(bl.Xgrid) points=length(bl.Y)
        return bl
    finally
        rm(tmp; recursive=true, force=true)
    end
end

"""
    load_ibl_solution(path) → NamedTuple | nothing

Reconstruct the IBL solution `(Xgrid, Y, u, v, w)` from an `ibl_out.mat` written
by the MATLAB driver. Internal helper for `run_ibl`; returns `nothing` if absent.
"""
function load_ibl_solution(path::AbstractString)
    isfile(path) || return nothing
    out = matread(path)
    return (Xgrid = vec(Float64.(out["Xgrid"])),
            Y     = vec(Float64.(out["Y"])),
            u     = Matrix{Float64}(out["u"]),
            v     = Matrix{Float64}(out["v"]),
            w     = Matrix{Float64}(out["w"]))
end
