"""
    make_tunnel_to_curved_plate(backend, root) → Dict{Symbol, Function}

Two-step pipeline: wind-tunnel simulation (TunnelCase), then mapped
curved-plate simulation (AirfoilLECase).
"""
function make_tunnel_to_curved_plate(backend::BackendType, root::AbstractString)
    tunnel_to_curved = joinpath(root, "PreProcessing", "Modules", "TunnelToCurvedPlateModule")

    tunnel_case   = joinpath(tunnel_to_curved, "TunnelCase")
    airfoil_case  = joinpath(tunnel_to_curved, "AirfoilLECase")
    generate_grid = joinpath(airfoil_case, "generateGrid.jl")
    map_script    = joinpath(tunnel_to_curved, "mapTunnelToAirfoilLE.jl")

    airfoil_data_dir = joinpath(root, "PreProcessing", "InputOutput", "AirfoilGeometryData")
    airfoil2stl      = joinpath(root, "PreProcessing", "Scripts", "SelfRunning", "airfoil2stl.sh")
    tri_surface_dir  = joinpath(tunnel_case, "constant", "triSurface")

    t = inp.TTCP.tunnel

    return Dict{Symbol, Function}(
        :clean => () -> begin
            foam_script(backend, tunnel_case, "clean")
            foam_script(backend, airfoil_case, "clean")
        end,

        :prep => () -> begin
            @info "Preparing airfoil STL..."
            mkpath(tri_surface_dir)

            dat_file = joinpath(airfoil_data_dir, t.airfoilFile)
            base     = replace(t.airfoilFile, ".dat" => "")
            raw_stl      = joinpath(tri_surface_dir, "$(base).stl")
            oriented_stl = joinpath(tri_surface_dir, "$(base)_oriented.stl")
            final_stl    = joinpath(tri_surface_dir, "airfoil.stl")

            # Tunnel z-extent [m] (hardcoded: z0=-0.01, z1=0.01 → 0.02 m)
            z_span = 0.02

            # Step 1: .dat → raw STL
            @info "  airfoil2stl: $(t.airfoilFile) → $(base).stl"
            run(`bash $airfoil2stl $dat_file $raw_stl $z_span`)

            # Step 2: Orient normals
            @info "  surfaceOrient → $(base)_oriented.stl"
            outside_pt = "($( t.tunnelLength / 2 ) 0 0)"
            foam_exec(backend, tunnel_case,
                      "surfaceOrient constant/triSurface/$(base).stl '$(outside_pt)' constant/triSurface/$(base)_oriented.stl")

            # Step 3: Check surface
            @info "  surfaceCheck..."
            foam_exec(backend, tunnel_case,
                      "surfaceCheck constant/triSurface/$(base)_oriented.stl")

            # Step 4: Scale, translate, rotate → airfoil.stl
            chord = t.chord
            dx = t.xCenter - 0.5 * chord
            dy = t.yCenter
            xo = t.xCenter
            yo = t.yCenter
            aoa = t.alphaDeg

            @info "  surfaceTransformPoints: chord=$(chord), alpha=$(aoa)°, center=($(t.xCenter),$(t.yCenter))"
            foam_exec(backend, tunnel_case,
                      "surfaceTransformPoints" *
                      " -scale '($chord $chord 1)'" *
                      " -translate '($dx 0 $dy)'" *
                      " -origin '($xo 0 $yo)'" *
                      " -rollPitchYaw '(0 0 $(-aoa))'" *
                      " constant/triSurface/$(base)_oriented.stl" *
                      " constant/triSurface/airfoil.stl")

            # Step 5: Extract feature edges
            @info "  surfaceFeatureExtract..."
            foam_exec(backend, tunnel_case, "surfaceFeatureExtract")

            # Clean up intermediate STLs
            rm(raw_stl;      force=true)
            rm(oriented_stl; force=true)

            @info "Airfoil STL ready."
        end,

        :meshTunnel => () -> begin
            write_tunnel_input_param(tunnel_case)

            @info "Meshing TunnelCase (blockMesh)..."
            foam_exec(backend, tunnel_case, "blockMesh")
            @info "Meshing TunnelCase (snappyHexMesh)..."
            foam_exec(backend, tunnel_case, "snappyHexMesh")
        end,

        :meshAirfoil => () -> begin
            write_airfoil_le_input_param(airfoil_case)

            @info "Generating AirfoilLECase grid..."
            run_julia_subprocess(generate_grid; dir=airfoil_case)
            @info "Meshing AirfoilLECase..."
            foam_exec(backend, airfoil_case, "blockMesh")
        end,

        :runTunnel => () -> begin
            write_tunnel_input_param(tunnel_case)

            @info "Solving TunnelCase..."
            foam_script(backend, tunnel_case, "run")
        end,

        :map => () -> begin
            @info "Mapping tunnel → airfoil..."
            run_julia_subprocess(map_script; dir=tunnel_to_curved)
        end,

        :runAirfoil => () -> begin
            write_airfoil_le_input_param(airfoil_case)

            @info "Solving AirfoilLECase..."
            foam_script(backend, airfoil_case, "run")
        end,

        :postAirfoil => () -> begin
            @info "Post-processing AirfoilLECase..."
            foam_script(backend, airfoil_case, "runPostProcess")
        end,

        :vizTunnel => () -> begin
            plotting_dir = joinpath(root, "PreProcessing", "InputOutput", "Plotting",
                                       "TunnelToCurvedPlate")
            return plot_residuals(tunnel_case; savedir=plotting_dir, label="TunnelCase")
        end,

        :vizAirfoil => () -> begin
            plotting_dir = joinpath(root, "PreProcessing", "InputOutput", "Plotting",
                                       "TunnelToCurvedPlate")

            res_airfoil = plot_residuals(airfoil_case; savedir=plotting_dir, label="AirfoilLECase")
            fields = plot_fields(airfoil_case; savedir=plotting_dir)
            bl = plot_profiles(airfoil_case; savedir=plotting_dir)

            return (res_airfoil=res_airfoil, fields=fields, bl=bl)
        end,

        :_order => () -> [:clean, :prep, :meshTunnel, :runTunnel,
                          :meshAirfoil, :map, :runAirfoil, :postAirfoil,
                          :vizTunnel, :vizAirfoil],
    )
end

MODULE_REGISTRY["TunnelToCurvedPlate"] = make_tunnel_to_curved_plate
