using Printf

@enum BackendType DOCKER NATIVE

"""
    detect_backend() → BackendType

macOS → DOCKER (OpenFOAM runs inside Docker container),
Linux → NATIVE (OpenFOAM installed natively).
"""
function detect_backend()
    Sys.isapple() ? DOCKER : NATIVE
end

# OpenFOAM bashrc path for NATIVE backend (Linux)
const OPENFOAM_BASHRC = let
    # Check if already sourced
    if haskey(ENV, "WM_PROJECT_DIR")
        joinpath(ENV["WM_PROJECT_DIR"], "etc", "bashrc")
    else
        # Search common install locations for any OpenFOAM version
        found = ""
        for base in ["/usr/lib/openfoam", "/opt", "/opt/OpenFOAM",
                     joinpath(homedir(), "OpenFOAM")]
            isdir(base) || continue
            for d in sort(readdir(base), rev=true)  # newest version first
                candidate = joinpath(base, d, "etc", "bashrc")
                if isfile(candidate)
                    found = candidate
                    break
                end
            end
            !isempty(found) && break
        end
        if isempty(found)
            @warn "OpenFOAM bashrc not found — source OpenFOAM before running, or set WM_PROJECT_DIR"
            "openfoam_not_found"
        else
            @info "Found OpenFOAM: $found"
            found
        end
    end
end

const DOCKER_CONTAINER = "openfoam2412"
const DOCKER_MOUNT_HOST = joinpath(homedir(), "OpenFOAM", "cases")
const DOCKER_MOUNT_CONTAINER = "/home/openfoam"

"""
    docker_case_path(host_path) → String

Translate a host path like `~/OpenFOAM/cases/airPower/...`
to the container path `/home/openfoam/airPower/...`.
"""
function docker_case_path(host_path::AbstractString)
    abs = abspath(host_path)
    if !startswith(abs, DOCKER_MOUNT_HOST)
        error("Path $abs is not under $DOCKER_MOUNT_HOST — cannot map to Docker container")
    end
    rel = abs[length(DOCKER_MOUNT_HOST)+2:end]
    return joinpath(DOCKER_MOUNT_CONTAINER, rel)
end

"""
    foam_exec(backend, case_path, cmd; verbose=true) → Bool

Run a single OpenFOAM command in `case_path`.
Returns true if the command succeeded.
"""
function foam_exec(backend::BackendType, case_path::AbstractString, cmd::AbstractString; verbose::Bool=true)
    if backend == DOCKER
        container_path = docker_case_path(case_path)
        full_cmd = `docker exec $(DOCKER_CONTAINER) bash -c "source /openfoam/bash.rc && cd $(container_path) && $(cmd)"`
    else
        full_cmd = `bash -c "cd $(abspath(case_path)) && source $(OPENFOAM_BASHRC) && $(cmd)"`
    end
    verbose && @info "Running: $cmd" case=basename(case_path)
    proc = run(ignorestatus(full_cmd), wait=true)
    if proc.exitcode != 0
        @error "Command failed (exit $(proc.exitcode)): $cmd"
        return false
    end
    return true
end

"""
    foam_script(backend, case_path, script; verbose=true) → Bool

Run an existing shell script (e.g. `run`, `clean`, `runPostProcess`) inside the case.
"""
function foam_script(backend::BackendType, case_path::AbstractString, script::AbstractString,
                     args::AbstractString=""; verbose::Bool=true)
    script_path = joinpath(case_path, script)
    if !isfile(script_path)
        @warn "Script not found: $script_path"
        return false
    end
    cmd = isempty(args) ? "bash $script" : "bash $script $args"
    foam_exec(backend, case_path, cmd; verbose)
end

# ── Wall modulation (smooth bump/depression) ─────────────────────────
include(joinpath(@__DIR__, "wallModulation.jl"))

"""
    write_flat_plate_input_param(case_dir)

Generate `constant/inputParam` for the FlatPlateModule from the central
Julia inputs.  The file is auto-generated and should not be edited by hand.
"""
function write_flat_plate_input_param(case_dir::AbstractString)
    p = inp.DFP
    wm = merge(inp.wallModulation, p.wallModulation)

    # Compute wall vertex y-displacements [mm] at the 7 block boundaries
    L = p.domainLength
    xverts = [0.0, 3.0/13.0*L, 7.0/26.0*L, 9.0/26.0*L, 5.0/13.0*L, 11.0/13.0*L, L]

    # Base cell counts per block (must match blockMeshDict base values)
    nx_base = [144, 24, 48, 24, 280, 96]
    gx = p.gridXfactor
    gy = p.gridYfactor

    if wm.enabled && wm.mode == :single
        # Vertices: evaluate bump (returns 0 outside full extent naturally)
        yverts = [wall_bump(xv, wm) * 1000.0 for xv in xverts]  # mm
    else
        yverts = zeros(7)
    end
    nx_final = [round(Int, nx_base[i] * gx) for i in 1:6]

    path = joinpath(case_dir, "constant", "inputParam")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
// AUTO-GENERATED by airPower orchestrator (inputs.jl) — do not edit
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      inputParam;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Boundaries information
Upinlet
{
//Domain geometry [m]
 domainLength  $(p.domainLength);
 domainHeight  $(p.domainHeight);

//Inflow parameters
 Uinf        $(p.Uinf);
 Winf        $(p.Winf);
 xInlet      $(p.xInlet);

//Top-boundary parameters (coefficient of polinomyal, Casacuberta et al, 2022)
 pa4         $(p.pa4);
 pa3         $(p.pa3);
 pa2         $(p.pa2);
 pa1         $(p.pa1);
 pa0         $(p.pa0);

//Fluid properties
 freeStreamViscosity   $(p.freeStreamViscosity);

//Grid resolution factors
 gridXfactor   $(p.gridXfactor);
 gridYfactor   $(p.gridYfactor);

//Output settings
 outputFormat          $(p.outputFormat);
 wallExtrapolation     $(p.wallExtrapolation ? "true" : "false");
}

//Wall modulation — vertex y-displacements [mm] at each block boundary
wallMod
{
 enabled  $(wm.enabled ? "true" : "false");
 y0       $(yverts[1]);
 y1       $(yverts[2]);
 y2       $(yverts[3]);
 y3       $(yverts[4]);
 y4       $(yverts[5]);
 y5       $(yverts[6]);
 y6       $(yverts[7]);
 Nx1      $(nx_final[1]);
 Nx2      $(nx_final[2]);
 Nx3      $(nx_final[3]);
 Nx4      $(nx_final[4]);
 Nx5      $(nx_final[5]);
 Nx6      $(nx_final[6]);
 NyB      $(round(Int, 120 * gy));
 NyT      $(round(Int, 40 * gy));
}

// ************************************************************************* //""")
    end
    @info "Generated constant/inputParam" case=basename(case_dir)

    # ── Write wall polyLine edges (system/wallEdges) ──
    edges_path = joinpath(case_dir, "constant", "wallEdges")
    N_KNOTS = 2000  # knots per edge segment (only affects blockMesh, not solver)

    open(edges_path, "w") do io
        println(io, "// AUTO-GENERATED wall modulation edges — do not edit")
        if wm.enabled && wm.mode == :single
            # Full bump extent (with blend) for polyLine edge overlap check
            if wm.shape == :sigmoidal
                bump_xStart = wm.xStart
                bump_xEnd   = wm.xEnd
            elseif wm.shape == :esn
                bump_xStart, bump_xEnd = _esn_geometry_full(wm)
            end

            # Wall edges: vertex pairs along the plate (z0 and z1 planes)
            # z0 plane: 0→1, 1→2, ..., 5→6  (vertices 0–6)
            # z1 plane: 21→22, ..., 26→27    (vertices 21–27)
            z0_mm = -1.0
            z1_mm =  1.0

            for (edge_i, (vi, vj)) in enumerate(zip(0:5, 1:6))
                xa = xverts[edge_i]
                xb = xverts[edge_i + 1]

                # Check if this edge overlaps with the bump
                if xb <= bump_xStart || xa >= bump_xEnd
                    continue  # entirely outside bump — straight edge is fine
                end

                # Generate knots: dense uniform within bump, sparse flat outside
                bump_lo = max(xa, bump_xStart)
                bump_hi = min(xb, bump_xEnd)

                # Uniform knots: y=0 outside bump, bump shape inside
                knots_x = collect(range(xa, xb, length=N_KNOTS + 2)[2:end-1])
                knots_y = [bump_xStart < xk < bump_xEnd ? wall_bump(xk, wm) * 1000.0 : 0.0
                           for xk in knots_x]
                knots_x_mm = knots_x .* 1000.0
                n_total = length(knots_x)

                # z0 plane
                println(io, "    polyLine $vi $(vi+1)")
                println(io, "    (")
                for k in 1:n_total
                    @Printf.printf(io, "        (%.10g %.10g %.10g)\n",
                                   knots_x_mm[k], knots_y[k], z0_mm)
                end
                println(io, "    )")

                # z1 plane (same x,y knots, different z)
                vi_z1 = vi + 21
                vj_z1 = vj + 21
                println(io, "    polyLine $vi_z1 $vj_z1")
                println(io, "    (")
                for k in 1:n_total
                    @Printf.printf(io, "        (%.10g %.10g %.10g)\n",
                                   knots_x_mm[k], knots_y[k], z1_mm)
                end
                println(io, "    )")
            end
        end
    end
    @info "Generated system/wallEdges" case=basename(case_dir)

    # ── Write refinement dicts for bump region (system/topoSetDict, system/refineMeshDict) ──
    if wm.enabled && wm.mode == :single
        # Refinement uses yTol extent (not the blend zones)
        if wm.shape == :sigmoidal
            bump_xs = wm.xStart
            bump_xe = wm.xEnd
        elseif wm.shape == :esn
            bump_xs, bump_xe = _esn_geometry(wm)
        end
        bumpL = bump_xe - bump_xs
        H = p.domainHeight

        # Three refinement boxes:
        #   box1 (outer, 2×): margin = bumpL/2 on each side
        #   box2 (middle, 4×): margin = bumpL/4 on each side
        #   box3 (inner, 8×): exactly the bump region
        x1_lo = max(0.0, bump_xs - bumpL / 2)
        x1_hi = min(L, bump_xe + bumpL / 2)
        x2_lo = max(0.0, bump_xs - bumpL / 4)
        x2_hi = min(L, bump_xe + bumpL / 4)
        x3_lo = bump_xs
        x3_hi = bump_xe

        # topoSetDict — defines three cellSets
        topo_path = joinpath(case_dir, "system", "topoSetDict")
        open(topo_path, "w") do io
            println(io, """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}

actions
(
    {
        name    refineBox1;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     ($x1_lo -1 -1) ($x1_hi $(H+1) 1);
    }
    {
        name    refineBox2;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     ($x2_lo -1 -1) ($x2_hi $(H+1) 1);
    }
    {
        name    refineBox3;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     ($x3_lo -1 -1) ($x3_hi $(H+1) 1);
    }
);""")
        end

        # refineMeshDict — refine in x only
        refine_path = joinpath(case_dir, "system", "refineMeshDict")
        open(refine_path, "w") do io
            println(io, """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      refineMeshDict;
}

set             refineBox1;

coordinateSystem global;

globalCoeffs
{
    tan1    (1 0 0);
    tan2    (0 1 0);
}

directions      (tan1);

useHexTopology  yes;

geometricCut    no;

writeMesh       no;""")
        end

        @info "Generated refinement dicts" case=basename(case_dir)
    end
end

const FOAM_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
// AUTO-GENERATED by airPower orchestrator (inputs.jl) — do not edit"""

"""
    write_decompose_par_dict(case_dir, nProcs)

Generate `system/decomposeParDict` with the given number of processors.
"""
function write_decompose_par_dict(case_dir::AbstractString, nProcs::Int)
    path = joinpath(case_dir, "system", "decomposeParDict")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, """$(FOAM_HEADER)
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains $(nProcs);

method          simple;

simpleCoeffs
{
    n               ( $(nProcs) 1 1 );
    delta           0.001;
}

// ************************************************************************* //""")
    end
end

"""
    write_tunnel_input_param(case_dir)

Generate `constant/inputParam` for TunnelCase from the central Julia inputs.
"""
function write_tunnel_input_param(case_dir::AbstractString)
    f = inp.TTCP.flow
    t = inp.TTCP.tunnel
    mm(x) = x * 1000.0   # m → mm
    path = joinpath(case_dir, "constant", "inputParam")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, """$(FOAM_HEADER)
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      inputParam;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

flowParam
{
//1.- Tunnel Geometrical Parameters

 //Basic geometry

 tunnelInletHalfHeight     $(mm(t.tunnelInletHalfHeight));        // [mm] (tunnel half height)

 tunnelOutletHalfHeight    $(mm(t.tunnelOutletHalfHeight));        // [mm]

 tunnelLength       $(mm(t.tunnelLength));         //[mm]

//2.- Grid Discretization Parameters

 //Streamwise partitions
 Nx $(t.Nx);                             // [-]

 //Wall-normal partitions
 Ny $(t.Ny);                              // [-]

//3.- Flow Parameters (m as dimensional units)

 freeStreamVelocityStreamwise    $(f.freeStreamVelocityStreamwise);       //[ms^1]

 freeStreamVelocitySpanwise     $(f.freeStreamVelocitySpanwise);        //[ms^1]

 freeStreamViscosity   $(f.freeStreamViscosity);    //[m^2s^1]

//4.- Turbulence Parameters

 turbulenceIntensity   $(t.turbulenceIntensity);        //[-]
 turbLengthScale       $(t.turbLengthScale);         //[m]

//5.- BSpline Control Points [mm]

 xcp1 $(mm(t.xcp1));
 ycp1 $(mm(t.ycp1));

 xcp2 $(mm(t.xcp2));
 ycp2 $(mm(t.ycp2));

 xcp3 $(mm(t.xcp3));
 ycp3 $(mm(t.ycp3));
}

// ************************************************************************* //""")
    end
    write_decompose_par_dict(case_dir, inp.TTCP.nProcs)
    @info "Generated constant/inputParam" case=basename(case_dir)
end

"""
    write_airfoil_le_input_param(case_dir)

Generate `constant/inputParam` for AirfoilLECase from the central Julia inputs.
(Grid/domain parameters in `system/inputDomain` are no longer needed —
 `generateGrid.jl` reads `inputs.jl` directly.)
"""
function write_airfoil_le_input_param(case_dir::AbstractString)
    f = inp.TTCP.flow
    t = inp.TTCP.tunnel
    a = inp.TTCP.airfoilLE
    mm(x) = x * 1000.0   # m → mm

    # Load upper-surface airfoil data from the .dat file
    airfoil_dat = joinpath(dirname(dirname(case_dir)),
                           "..", "InputOutput", "AirfoilGeometryData", t.airfoilFile)
    xi_up, eta_up = load_airfoil_upper(airfoil_dat)

    # ── constant/inputParam ──────────────────────────────────────────────
    path_ip = joinpath(case_dir, "constant", "inputParam")
    mkpath(dirname(path_ip))
    open(path_ip, "w") do io
        println(io, """$(FOAM_HEADER)
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      inputParam;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Boundaries information
Upinlet
{
//Fluid properties
 freeStreamViscosity            $(f.freeStreamViscosity);              //[m^2s^1]
 freeStreamVelocityStreamwise   $(f.freeStreamVelocityStreamwise);    //[ms^1]
 freeStreamVelocitySpanwise     $(f.freeStreamVelocitySpanwise);      //[ms^1]

//Output settings
 outputFormat          $(a.outputFormat);      // csv | binary
 wallExtrapolation     $(a.wallExtrapolation ? "true" : "false");

//Export region settings
 exportMode            $(a.exportMode);  // full | partial
 xiInlet               $(a.xiInlet);     // chord fraction for inlet boundary (x/c)
 xiOutlet              $(a.xiOutlet);     // chord fraction for outlet boundary (x/c)
 exportHeight          $(mm(a.exportHeight));       // wall-normal distance from surface [mm]

//Airfoil geometry
 chord                 $(mm(t.chord));      // [mm]
 alphaDeg              $(t.alphaDeg);     // [deg]
 xCenter               $(mm(t.xCenter));     // [mm]
 yCenter               $(mm(t.yCenter));     // [mm]

//Upper-surface airfoil coordinates (xi/c, eta/c), LE to TE
 airfoilXi             $(length(xi_up))($(join(string.(xi_up), " ")));
 airfoilEta            $(length(eta_up))($(join(string.(eta_up), " ")));
}

//Boundary-layer integral metrics (consumed by blMetrics.jl post-processor)
blMetrics
{
 method                $(a.blMetrics.method);   // vorticityIntegralTrapezoidal | vorticityIntegralMidpoint | maxProfile | fixedHeight | pressureBernoulli
}

// ************************************************************************* //""")
    end
    write_decompose_par_dict(case_dir, inp.TTCP.nProcs)
    @info "Generated constant/inputParam" case=basename(case_dir)
end

"""
    load_airfoil_upper(filepath) → (xi, eta)

Read an airfoil .dat file (tab/space-separated x/c, y/c columns, full
contour TE→LE→TE) and return the upper-surface points sorted LE→TE.
"""
function load_airfoil_upper(filepath::AbstractString)
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

    mask = eta_all .>= 0.0
    xi_up  = xi_all[mask]
    eta_up = eta_all[mask]

    perm = sortperm(xi_up)
    xi_up  = xi_up[perm]
    eta_up = eta_up[perm]

    keep = [true; [!(xi_up[i] == xi_up[i-1] && eta_up[i] == eta_up[i-1]) for i in 2:length(xi_up)]]
    return xi_up[keep], eta_up[keep]
end

"""
    run_julia_subprocess(script_path; dir, args) → Bool

Run a Julia script as a subprocess (isolated from our module scope).
"""
function run_julia_subprocess(script_path::AbstractString;
                              dir::AbstractString=dirname(script_path),
                              args::Vector{String}=String[])
    abs_script = abspath(script_path)
    if !isfile(abs_script)
        @warn "Julia script not found: $abs_script"
        return false
    end
    @info "Running Julia subprocess: $(basename(abs_script))"
    cmd = `julia $abs_script $args`
    proc = run(Cmd(cmd; dir), wait=true)
    return proc.exitcode == 0
end
