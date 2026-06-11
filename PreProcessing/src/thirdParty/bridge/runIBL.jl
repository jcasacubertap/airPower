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

# blockMesh streamwise base cell counts per block (must match the nx_base in
# backend.jl write_flat_plate_input_param). Total × gridXfactor is the
# streamwise cell count of the DFP mesh.
const _DFP_NX_BASE = (144, 24, 48, 24, 280, 96)   # sum = 616

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
    run_ibl_solver(case_path; savedir) → NamedTuple | nothing

Run the spectral integral boundary-layer solver (solver_IBL_spectral, via
MATLAB) for the DirectFlatPlate case and return its velocity fields for overlay
on the DFP plots. Returns `nothing` if the solver cannot be run (e.g. MATLAB
unavailable) so that visualization can proceed without the overlay.

The spectral solver seeds u, v AND w at the inlet from the Falkner-Skan-Cooke
similarity solution, so the crossflow is imposed at x = 0 and matches the DFP
at every station (the Thomas-algorithm solver_IBL grows w from a uniform inflow
and does not).

Solver inputs are derived from `inp.DFP` and `inp.VAL`:
  S  = xInlet,  L = xInlet + domainLength,  H = domainHeight,
  nu = freeStreamViscosity,  We = Winf,
  nx  = blockMesh streamwise cell count (sum(_DFP_NX_BASE)·gridXfactor),
  ny  = inp.VAL.blNcheb (number of wall-normal Chebyshev nodes),
  y_i = inp.VAL.blYi (Chebyshev node median — half the nodes below it),
  Ue  = analytic edge velocity `Uinf·P(log x)` — the SAME distribution the DFP
        imposes (see `_ue_poly`), so the inlet FSC beta coincides with the DFP.

`case_path` is retained for signature compatibility / future field checks; the
solver inputs do not depend on parsing the DFP field output.

Returned NamedTuple fields:
  Xgrid  streamwise coordinate in the DFP frame (0 -> domainLength)
  Y      wall-normal Chebyshev coordinate [m] (H -> 0)
  u,v,w  ny x nxout IBL velocity components (row 1 = freestream, row ny = wall)
"""
function run_ibl_solver(case_path::AbstractString;
                        savedir::AbstractString=case_path)
    matlab = _find_matlab()
    if matlab === nothing
        @warn "MATLAB executable not found — skipping IBL overlay (set inp.VAL.valBL=false to silence)."
        return nothing
    end

    p = inp.DFP
    S = p.xInlet
    L = p.xInlet + p.domainLength
    H = p.domainHeight
    nu = p.freeStreamViscosity
    We = p.Winf

    # Streamwise resolution from the blockMesh. Wall-normal grid is the
    # spectral solver's Chebyshev grid: blNcheb nodes mapped over [0, H] with
    # median node at y_i = blYi (half the nodes below it — cluster in the BL).
    nx  = sum(_DFP_NX_BASE) * Int(p.gridXfactor)
    ny  = hasproperty(inp.VAL, :blNcheb) ? Int(inp.VAL.blNcheb) : 150
    y_i = hasproperty(inp.VAL, :blYi)    ? Float64(inp.VAL.blYi) : 0.002

    # Edge velocity: the SAME analytic distribution the DFP imposes (log-
    # pressure polynomial), sampled on the uniform streamwise grid x_OF in
    # [0, domainLength] → physical arc-length x = S .. L. Because Ue(x) (and its
    # inlet gradient) matches the DFP, the solver derives the same FSC beta and
    # the imposed inlet profile coincides with the DFP — no solver change.
    xg = range(0.0, p.domainLength, length=nx)
    Ue = [_ue_poly(S + xq, p) for xq in xg]

    @info "IBL spectral grid & edge velocity" nx=nx ny_cheb=ny y_i=y_i H=H Ue_inlet=round(Ue[1], digits=4) Ue_outlet=round(Ue[end], digits=4)

    driver_dir = @__DIR__                          # PreProcessing/src/thirdParty/bridge
    mkpath(savedir)
    infile  = joinpath(savedir, "bl_input.mat")
    outfile = joinpath(savedir, "bl_output.mat")
    isfile(outfile) && rm(outfile; force=true)

    matwrite(infile, Dict(
        "nx"  => Float64(nx), "ny" => Float64(ny),
        "S"   => Float64(S),  "L"  => Float64(L), "H" => Float64(H),
        "y_i" => Float64(y_i),
        "nu"  => Float64(nu), "We" => Float64(We),
        "Ue"  => Vector{Float64}(Ue),
    ))

    cmd = Cmd(`$matlab -batch "run_IBL_DFP('$infile','$outfile')"`; dir=driver_dir)
    @info "Invoking MATLAB IBL driver..." matlab=matlab
    try
        run(cmd)
    catch e
        @warn "MATLAB IBL driver failed — skipping overlay." exception=e
        return nothing
    end

    if !isfile(outfile)
        @warn "IBL driver produced no output ($outfile) — skipping overlay."
        return nothing
    end

    out = matread(outfile)
    Xgrid = vec(Float64.(out["Xgrid"]))
    Y     = vec(Float64.(out["Y"]))
    u = Matrix{Float64}(out["u"])
    v = Matrix{Float64}(out["v"])
    w = Matrix{Float64}(out["w"])
    @info "IBL solution loaded" stations=length(Xgrid) points=length(Y)
    return (Xgrid=Xgrid, Y=Y, u=u, v=v, w=w)
end
