#!/usr/bin/env julia
#=  ═══════════════════════════════════════════════════════════════════════════
    mapTunnelToAirfoilLE.jl

    Maps the converged TunnelCase RANS solution onto AirfoilLECase boundaries.

    Workflow:
      1. Parse AirfoilLECase mesh to extract face centers for target patches
      2. Write an OpenFOAM sampleDict (cloud-type sets at those locations)
      3. Run postProcess on TunnelCase to sample U and p
      4. Write constant/boundaryData/ files for timeVaryingMappedFixedValue
      5. Update 0/U and 0/p boundary conditions

    Prerequisites:
      - AirfoilLECase mesh must already exist (run generateGrid.jl && blockMesh)
      - TunnelCase must be converged (time directory must exist)
      - OpenFOAM must be sourced in the environment
    ═══════════════════════════════════════════════════════════════════════════ =#

using Printf

# ═══════════════════════════════════════════════════════════════════════════════
#  Configuration
# ═══════════════════════════════════════════════════════════════════════════════

const SCRIPT_DIR   = @__DIR__
const TUNNEL_CASE  = joinpath(SCRIPT_DIR, "TunnelCase")
const AIRFOIL_CASE = joinpath(SCRIPT_DIR, "AirfoilLECase")

const TUNNEL_TIME = "163"   # converged time step
const Z_PROBE     = 0.005   # z [m] for sampling (TunnelCase cell centre, ±5 mm)

# ═══════════════════════════════════════════════════════════════════════════════
#  OpenFOAM mesh parsing
# ═══════════════════════════════════════════════════════════════════════════════

"""
    skip_foam_header(io) -> nothing

Advance the stream past the FoamFile header and the `// * * ...` separator line.
"""
function skip_foam_header(io::IO)
    for line in eachline(io)
        if startswith(line, "// *")
            return
        end
    end
end

"""
    parse_points(filepath) -> Vector{NTuple{3,Float64}}

Read an OpenFOAM `points` file (vectorField).
"""
function parse_points(filepath::String)
    open(filepath) do io
        skip_foam_header(io)

        # Read until we find the count line (a bare integer)
        n = 0
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            m = match(r"^(\d+)$", s)
            if m !== nothing
                n = parse(Int, m.captures[1])
                break
            end
        end
        @assert n > 0 "Could not find point count in $filepath"

        # Skip the opening '('
        while true
            s = strip(readline(io))
            s == "(" && break
        end

        points = Vector{NTuple{3,Float64}}(undef, n)
        for i in 1:n
            line = readline(io)
            # Format: (x y z)
            m = match(r"\(\s*([^ ]+)\s+([^ ]+)\s+([^ )]+)\s*\)", line)
            points[i] = (parse(Float64, m.captures[1]),
                         parse(Float64, m.captures[2]),
                         parse(Float64, m.captures[3]))
        end
        return points
    end
end

"""
    parse_boundary(filepath) -> Dict{String, NamedTuple{(:nFaces,:startFace), ...}}

Read an OpenFOAM `boundary` file and return patch info.
"""
function parse_boundary(filepath::String)
    text = read(filepath, String)
    patches = Dict{String, NamedTuple{(:nFaces, :startFace), Tuple{Int,Int}}}()

    # Find each patch block: name { ... nFaces N; startFace S; ... }
    for m in eachmatch(r"(\w+)\s*\{([^}]+)\}", text)
        name = m.captures[1]
        block = m.captures[2]

        mf = match(r"nFaces\s+(\d+)", block)
        ms = match(r"startFace\s+(\d+)", block)
        if mf !== nothing && ms !== nothing
            patches[name] = (nFaces    = parse(Int, mf.captures[1]),
                             startFace = parse(Int, ms.captures[1]))
        end
    end
    return patches
end

"""
    parse_faces_range(filepath, start_face, n_faces) -> Vector{Vector{Int}}

Read `n_faces` face entries starting at index `start_face` (0-based face index).
Returns vertex indices (0-based, matching OpenFOAM convention).
"""
function parse_faces_range(filepath::String, start_face::Int, n_faces::Int)
    faces = Vector{Vector{Int}}(undef, n_faces)

    open(filepath) do io
        skip_foam_header(io)

        # Find the count
        while true
            s = strip(readline(io))
            isempty(s) && continue
            m = match(r"^(\d+)$", s)
            if m !== nothing
                break
            end
        end

        # Skip the opening '('
        while true
            s = strip(readline(io))
            s == "(" && break
        end

        # Skip to start_face
        for _ in 1:start_face
            readline(io)
        end

        # Parse n_faces lines
        for i in 1:n_faces
            line = strip(readline(io))
            # Format: N(v1 v2 ... vN)
            m = match(r"\d+\(([^)]+)\)", line)
            faces[i] = parse.(Int, split(m.captures[1]))
        end
    end
    return faces
end

"""
    compute_face_centers(faces, points) -> Vector{NTuple{3,Float64}}

Compute face centres as the average of their vertex positions.
"""
function compute_face_centers(faces::Vector{Vector{Int}},
                              points::Vector{NTuple{3,Float64}})
    centers = Vector{NTuple{3,Float64}}(undef, length(faces))
    for (i, face) in enumerate(faces)
        cx, cy, cz = 0.0, 0.0, 0.0
        for vi in face
            p = points[vi + 1]  # 0-based → 1-based
            cx += p[1]; cy += p[2]; cz += p[3]
        end
        nv = length(face)
        centers[i] = (cx / nv, cy / nv, cz / nv)
    end
    return centers
end

# ═══════════════════════════════════════════════════════════════════════════════
#  sampleDict generation
# ═══════════════════════════════════════════════════════════════════════════════

function write_sample_dict(patch_centers::Dict{String, Vector{NTuple{3,Float64}}})
    dict_path = joinpath(TUNNEL_CASE, "system", "sampleDict")

    open(dict_path, "w") do io
        print(io, """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type            sets;
libs            (sampling);

interpolationScheme cellPoint;
setFormat       raw;

fields          (U p);

sets
{
""")

        for (patch_name, centers) in patch_centers
            set_name = patch_name * "Probes"
            println(io, "    $set_name")
            println(io, "    {")
            println(io, "        type    cloud;")
            println(io, "        axis    xyz;")
            println(io, "        points")
            println(io, "        (")
            for c in centers
                @printf(io, "            (%.10g %.10g %.10g)\n",
                        c[1], c[2], Z_PROBE)
            end
            println(io, "        );")
            println(io, "    }")
        end

        println(io, "}")
        println(io)
        println(io, "// ************************************************************************* //")
    end

    println("  Written: $dict_path")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Run OpenFOAM sampling
# ═══════════════════════════════════════════════════════════════════════════════

"""
    docker_case_path(host_path)

Convert a host path under ~/OpenFOAM/cases/ to the Docker mount point
/home/openfoam/ used by the opencfd/openfoam-default container.
"""
function docker_case_path(host_path::String)
    home = ENV["HOME"]
    prefix = joinpath(home, "OpenFOAM", "cases")
    if startswith(host_path, prefix)
        return "/home/openfoam" * host_path[length(prefix)+1:end]
    end
    error("Path $host_path is not under $prefix")
end

function run_sampling()
    container_case = docker_case_path(TUNNEL_CASE)
    foam_cmd = "postProcess -case $container_case -func sampleDict -time $TUNNEL_TIME"

    println("  Running inside Docker container openfoam2412:")
    println("    $foam_cmd")

    run(`docker exec openfoam2412 bash -ic $foam_cmd`)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Parse sampled data (.xy files from the raw set format)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    parse_sampled_xy(filepath)

Parse a combined raw `.xy` file produced by the sample utility.
Columns: x y z p Ux Uy Uz   (all sampled fields concatenated).
Returns (p_values, U_values).
"""
function parse_sampled_xy(filepath::String)
    lines = filter(!isempty, strip.(readlines(filepath)))
    n = length(lines)

    p_values = Vector{Float64}(undef, n)
    U_values = Vector{NTuple{3,Float64}}(undef, n)

    for (i, line) in enumerate(lines)
        cols = split(line)
        p_values[i] = parse(Float64, cols[4])
        U_values[i] = (parse(Float64, cols[5]),
                       parse(Float64, cols[6]),
                       parse(Float64, cols[7]))
    end
    return p_values, U_values
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Write boundaryData for timeVaryingMappedFixedValue
# ═══════════════════════════════════════════════════════════════════════════════

function write_boundary_data(patch_name::String,
                             face_centers::Vector{NTuple{3,Float64}},
                             field_name::String,
                             values,
                             is_vector::Bool)
    bd_dir   = joinpath(AIRFOIL_CASE, "constant", "boundaryData", patch_name)
    time_dir = joinpath(bd_dir, "0")
    mkpath(time_dir)

    n = length(face_centers)

    # ── points file ──
    points_path = joinpath(bd_dir, "points")
    open(points_path, "w") do io
        println(io, n)
        println(io, "(")
        for c in face_centers
            @printf(io, "    (%.10g %.10g %.10g)\n", c[1], c[2], c[3])
        end
        println(io, ")")
    end
    println("  Written: $points_path  ($n points)")

    # ── field file ──
    field_path = joinpath(time_dir, field_name)
    open(field_path, "w") do io
        println(io, n)
        println(io, "(")
        if is_vector
            for v in values
                @printf(io, "    (%.10g %.10g %.10g)\n", v[1], v[2], v[3])
            end
        else
            for v in values
                @printf(io, "    %.10g\n", v)
            end
        end
        println(io, ")")
    end
    println("  Written: $field_path")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Update AirfoilLECase boundary conditions
# ═══════════════════════════════════════════════════════════════════════════════

function update_boundary_conditions()
    # ── 0/U ──
    u_path = joinpath(AIRFOIL_CASE, "0", "U")
    u_text = read(u_path, String)

    # Remove #include and Uinf lines if present
    u_text = replace(u_text, r"#include\s+\"[^\"]*inputParam\"\s*\n" => "")
    u_text = replace(u_text, r"Uinf\s+[^;]+;\s*\n" => "")

    # Set internalField to uniform (0 0 0) — potentialFoam will overwrite
    u_text = replace(u_text,
        r"internalField\s+[^;]+;" =>
        "internalField   uniform (0 0 0);")

    # farfield: fixedValue → timeVaryingMappedFixedValue
    u_text = replace(u_text,
        r"farfield\s*\{[^}]*\}" =>
        """farfield
    {
        type            timeVaryingMappedFixedValue;
        offset          (0 0 0);
        setAverage      false;
        value           uniform (0 0 0);
    }""")

    # suctionOutlet: zeroGradient → inletOutlet
    u_text = replace(u_text,
        r"suctionOutlet\s*\{[^}]*\}" =>
        """suctionOutlet
    {
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           uniform (0 0 0);
    }""")

    # pressureOutlet: zeroGradient → inletOutlet
    u_text = replace(u_text,
        r"pressureOutlet\s*\{[^}]*\}" =>
        """pressureOutlet
    {
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           uniform (0 0 0);
    }""")

    write(u_path, u_text)
    println("  Updated: $u_path")

    # ── 0/p ──
    p_path = joinpath(AIRFOIL_CASE, "0", "p")
    p_text = read(p_path, String)

    # suctionOutlet: fixedValue → timeVaryingMappedFixedValue
    p_text = replace(p_text,
        r"suctionOutlet\s*\{[^}]*\}" =>
        """suctionOutlet
    {
        type            timeVaryingMappedFixedValue;
        offset          0;
        setAverage      false;
        value           uniform 0;
    }""")

    # pressureOutlet: fixedValue → timeVaryingMappedFixedValue
    p_text = replace(p_text,
        r"pressureOutlet\s*\{[^}]*\}" =>
        """pressureOutlet
    {
        type            timeVaryingMappedFixedValue;
        offset          0;
        setAverage      false;
        value           uniform 0;
    }""")

    write(p_path, p_text)
    println("  Updated: $p_path")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

function main()
    println("╔══════════════════════════════════════════════════════════════╗")
    println("║  Map TunnelCase → AirfoilLECase boundary data              ║")
    println("╚══════════════════════════════════════════════════════════════╝")
    println()

    # Patch → which fields to write as boundaryData
    patch_fields = Dict(
        "farfield"       => ["U"],
        "suctionOutlet"  => ["p"],
        "pressureOutlet" => ["p"],
    )

    # All patches to sample
    patch_names = collect(keys(patch_fields))

    # ── Step 1: Parse AirfoilLECase mesh ─────────────────────────────────────
    println("[Step 1] Parsing AirfoilLECase mesh...")

    points_file   = joinpath(AIRFOIL_CASE, "constant", "polyMesh", "points")
    faces_file    = joinpath(AIRFOIL_CASE, "constant", "polyMesh", "faces")
    boundary_file = joinpath(AIRFOIL_CASE, "constant", "polyMesh", "boundary")

    println("  Reading points...")
    points = parse_points(points_file)
    println("  Read $(length(points)) points")

    println("  Reading boundary...")
    boundary = parse_boundary(boundary_file)

    patch_centers = Dict{String, Vector{NTuple{3,Float64}}}()

    for patch_name in keys(patch_fields)
        info = boundary[patch_name]
        println("  Reading faces for $patch_name " *
                "($(info.nFaces) faces from startFace=$(info.startFace))...")
        faces   = parse_faces_range(faces_file, info.startFace, info.nFaces)
        centers = compute_face_centers(faces, points)
        patch_centers[patch_name] = centers
        println("    → $(length(centers)) face centres extracted")
    end
    println()

    # ── Step 2: Write sampleDict ─────────────────────────────────────────────
    println("[Step 2] Writing sampleDict for TunnelCase...")
    write_sample_dict(patch_centers)
    println()

    # ── Step 3: Run sampling on TunnelCase ───────────────────────────────────
    println("[Step 3] Sampling TunnelCase at t=$TUNNEL_TIME...")
    run_sampling()
    println()

    # ── Step 4: Parse sampled data and write boundaryData ────────────────────
    println("[Step 4] Writing boundaryData...")

    pp_dir = joinpath(TUNNEL_CASE, "postProcessing", "sampleDict", TUNNEL_TIME)

    for (patch_name, fields) in patch_fields
        set_name = patch_name * "Probes"
        centers  = patch_centers[patch_name]

        # Combined file: <setName>_p_U.xy  (columns: x y z p Ux Uy Uz)
        xy_file = joinpath(pp_dir, "$(set_name)_p_U.xy")
        println("  Parsing $xy_file...")
        p_values, U_values = parse_sampled_xy(xy_file)

        @assert length(p_values) == length(centers) (
            "Mismatch: $(length(p_values)) sampled values vs " *
            "$(length(centers)) face centres for $patch_name")

        for field in fields
            if field == "U"
                write_boundary_data(patch_name, centers, "U", U_values, true)
            else
                write_boundary_data(patch_name, centers, "p", p_values, false)
            end
        end
    end
    println()

    # ── Step 5: Update boundary conditions ───────────────────────────────────
    println("[Step 5] Updating AirfoilLECase boundary conditions...")
    update_boundary_conditions()
    println()

    println("Done!")
    println("  Next steps:")
    println("    1. Inspect constant/boundaryData/ files")
    println("    2. Run AirfoilLECase solver")
end

main()
