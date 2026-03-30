"""
    make_direct_flat_plate(backend, root) → Dict{Symbol, Function}

Direct flat-plate base-flow computation (single OpenFOAM case).
"""
function make_direct_flat_plate(backend::BackendType, root::AbstractString)
    case_dir = joinpath(root, "PreProcessing", "Modules",
                        "DirectFlatPlateModule")

    flow_data_dir = joinpath(root, "PreProcessing", "InputOutput", "AirfoilFlowData")
    plotting_dir  = joinpath(root, "PreProcessing", "InputOutput", "Plotting",
                             "DirectFlatPlate")

    return Dict{Symbol, Function}(
        :clean => () -> foam_script(backend, case_dir, "clean"),

        :prep => () -> begin
            v = inp.VAL.externalToScaling

            # Load flow data from .mat file
            mat_path = joinpath(flow_data_dir, v.flowDataFile)
            @info "Loading flow data" file=v.flowDataFile
            BL = matread(mat_path)["BL"]

            # Run external-to-scaling
            result = external_to_scaling(
                x            = vec(BL["x"]),
                S            = vec(BL["S"]),
                expUe        = vec(BL["Ue"]),
                chord        = BL["c"],
                percent_crop = v.percentCrop,
                Nfit         = v.Nfit,
                fit_law      = v.fitLaw,
                nuinf        = BL["nu"],
                savedir      = plotting_dir,
            )

            @info "Scaling properties at (virtual inlet) airfoil chord %" delta0=result.delta0 uinf=result.uinf beta_FSK=result.beta_FSK arclength=(result.S_crop[1])
            @info "Polynomial coefficients (highest order first):"
            for (i, c) in enumerate(result.ue_pol_coeff)
                @info "  a$(length(result.ue_pol_coeff)-i) = $c"
            end
            return result
        end,

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
            res    = plot_residuals(case_dir; savedir=plotting_dir, label="DirectFlatPlate")
            fields = plot_fields(case_dir; savedir=plotting_dir)

            return (residuals=res, fields=fields)
        end,

        :_order => () -> [:clean, :prep, :mesh, :solve, :post, :viz],
    )
end

MODULE_REGISTRY["DirectFlatPlate"] = make_direct_flat_plate
