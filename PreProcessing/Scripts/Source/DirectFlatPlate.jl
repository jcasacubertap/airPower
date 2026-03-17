"""
    make_direct_flat_plate(backend, root) → Dict{Symbol, Function}

Direct flat-plate base-flow computation (single OpenFOAM case).
"""
function make_direct_flat_plate(backend::BackendType, root::AbstractString)
    case_dir = joinpath(root, "PreProcessing", "Modules",
                        "DirectFlatPlateModule")

    return Dict{Symbol, Function}(
        :clean => () -> foam_script(backend, case_dir, "clean"),
        :mesh  => () -> begin
            write_flat_plate_input_param(case_dir)
            foam_exec(backend, case_dir, "blockMesh")
        end,
        :solve => () -> begin
            write_flat_plate_input_param(case_dir)
            foam_script(backend, case_dir, "run", "$(inp.DFP.nProcs)")
        end,

        :post => () -> begin
            @info "Post-processing DirectFlatPlate..."
            foam_script(backend, case_dir, "runPostProcess", "$(inp.DFP.nProcs)")
        end,

        :viz => () -> begin
            plotting_dir = joinpath(root, "PreProcessing", "InputOutput", "Plotting",
                                       "DirectFlatPlate")
            res    = plot_residuals(case_dir; savedir=plotting_dir, label="DirectFlatPlate")
            fields = plot_fields(case_dir; savedir=plotting_dir)

            return (residuals=res, fields=fields)
        end,
    )
end

MODULE_REGISTRY["DirectFlatPlate"] = make_direct_flat_plate
