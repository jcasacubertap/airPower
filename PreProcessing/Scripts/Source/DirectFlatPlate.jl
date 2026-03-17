"""
    make_direct_flat_plate(backend, root) → Dict{Symbol, Function}

Direct flat-plate base-flow computation (single OpenFOAM case).
"""
function make_direct_flat_plate(backend::BackendType, root::AbstractString)
    case_dir = joinpath(root, "PreProcessing", "Modules", "BaseFlowGenerator",
                        "DirectFlatPlateModule")

    return Dict{Symbol, Function}(
        :clean => () -> foam_script(backend, case_dir, "clean"),
        :mesh  => () -> begin
            write_flat_plate_input_param(case_dir)
            foam_exec(backend, case_dir, "blockMesh")
        end,
        :solve => () -> begin
            write_flat_plate_input_param(case_dir)
            foam_script(backend, case_dir, "run")
        end,

        :post => () -> begin
            @info "Post-processing midPlane (serial)..."
            foam_exec(backend, case_dir, "simpleFoam -postProcess -latestTime")
        end,

        :viz => () -> begin
            plotting_dir = joinpath(root, "PreProcessing", "InputOutput", "Plotting",
                                       "BaseFlowGenerator", "DirectFlatPlate")
            res    = plot_residuals(case_dir; savedir=plotting_dir, label="DirectFlatPlate")
            fields = plot_fields(case_dir; savedir=plotting_dir)

            return (residuals=res, fields=fields)
        end,
    )
end

MODULE_REGISTRY["DirectFlatPlate"] = make_direct_flat_plate
