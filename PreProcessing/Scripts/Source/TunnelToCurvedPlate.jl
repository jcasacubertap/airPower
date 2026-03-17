"""
    make_tunnel_to_curved_plate(backend, root) → Dict{Symbol, Function}

Two-step pipeline: wind-tunnel simulation (TunnelCase), then mapped
curved-plate simulation (AirfoilLECase).
"""
function make_tunnel_to_curved_plate(backend::BackendType, root::AbstractString)
    base = joinpath(root, "PreProcessing", "Modules", "BaseFlowGenerator")
    tunnel_to_curved = joinpath(base, "TunnelToCurvedPlateModule")

    tunnel_case   = joinpath(tunnel_to_curved, "TunnelCase")
    airfoil_case  = joinpath(tunnel_to_curved, "AirfoilLECase")
    generate_grid = joinpath(airfoil_case, "generateGrid.jl")
    map_script    = joinpath(tunnel_to_curved, "mapTunnelToAirfoilLE.jl")

    return Dict{Symbol, Function}(
        :clean => () -> begin
            foam_script(backend, tunnel_case, "clean")
            foam_script(backend, airfoil_case, "clean")
        end,

        :mesh => () -> begin
            @info "Meshing TunnelCase..."
            foam_exec(backend, tunnel_case, "blockMesh")
            @info "Generating AirfoilLECase grid..."
            run_julia_subprocess(generate_grid; dir=airfoil_case)
            @info "Meshing AirfoilLECase..."
            foam_exec(backend, airfoil_case, "blockMesh")
        end,

        :solve => () -> begin
            @info "Solving TunnelCase..."
            foam_script(backend, tunnel_case, "run")
            @info "Mapping tunnel → airfoil..."
            run_julia_subprocess(map_script; dir=tunnel_to_curved)
            @info "Solving AirfoilLECase..."
            foam_script(backend, airfoil_case, "run")
        end,

        :post => () -> begin
            @info "Post-processing AirfoilLECase..."
            foam_script(backend, airfoil_case, "runPostProcess")
        end,

        :viz => () -> begin
            plotting_dir = joinpath(root, "PreProcessing", "InputOutput", "Plotting",
                                       "BaseFlowGenerator", "TunnelToCurvedPlate")

            # Residuals for both cases
            res_tunnel  = plot_residuals(tunnel_case;  savedir=plotting_dir, label="TunnelCase")
            res_airfoil = plot_residuals(airfoil_case; savedir=plotting_dir, label="AirfoilLECase")

            # Field contours from AirfoilLECase midPlane
            fields = plot_fields(airfoil_case; savedir=plotting_dir)

            # BL profiles
            bl = plot_profiles(airfoil_case; savedir=plotting_dir)

            return (res_tunnel=res_tunnel, res_airfoil=res_airfoil,
                    fields=fields, bl=bl)
        end,
    )
end

MODULE_REGISTRY["TunnelToCurvedPlate"] = make_tunnel_to_curved_plate
