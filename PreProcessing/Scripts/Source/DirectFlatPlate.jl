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
            @info "Polynomial coefficients Ue/Uinf (highest order first):"
            for (i, c) in enumerate(result.ue_pol_coeff)
                @info "  a$(length(result.ue_pol_coeff)-i) = $(c / result.uinf)"
            end
            return result
        end,

        :mesh  => () -> begin
            write_flat_plate_input_param(case_dir)
            @info "Meshing DirectFlatPlate..."
            if backend == DOCKER
                work = "/tmp/DFP_mesh"
                foam_exec(backend, case_dir,
                    "rm -rf $work && cp -r . $work && cd $work && blockMesh")
                dfp_docker = docker_case_path(case_dir)
                run(ignorestatus(`docker exec $(DOCKER_CONTAINER) bash -c
                    "rm -rf $dfp_docker/constant/polyMesh && cp -r $work/constant/polyMesh $dfp_docker/constant/polyMesh"`))
                foam_exec(backend, case_dir, "rm -rf $work")
            else
                foam_exec(backend, case_dir, "blockMesh")
            end
        end,

        :run => () -> begin
            write_flat_plate_input_param(case_dir)
            @info "Solving DirectFlatPlate..."
            if backend == DOCKER
                work = "/tmp/DFP_run"
                np = inp.DFP.nProcs
                foam_exec(backend, case_dir,
                    "rm -rf $work && cp -r . $work && cd $work" *
                    " && decomposePar" *
                    " && mpirun -np $np --oversubscribe simpleFoam -parallel" *
                    " && reconstructPar" *
                    " && for f in FSC_inletData.dat FSC_inletProfile.dat; do" *
                    "   [ -f processor0/\$f ] && cp processor0/\$f .; done")
                latest = read(`docker exec $(DOCKER_CONTAINER) bash -c
                    "ls -1d $work/[0-9]* 2>/dev/null | grep -v '/0\$' | sort -g | tail -1"`, String) |> strip
                if !isempty(latest)
                    run(`docker cp $(DOCKER_CONTAINER):$latest $(case_dir)/`)
                    @info "Copied time directory: $(basename(latest))"
                end
                for f in ["FSC_inletData.dat", "FSC_inletProfile.dat"]
                    run(ignorestatus(`docker cp $(DOCKER_CONTAINER):$work/$f $(case_dir)/`))
                end
                foam_exec(backend, case_dir, "rm -rf $work")
            else
                foam_script(backend, case_dir, "run", "$(inp.DFP.nProcs)")
            end
        end,

        :post => () -> begin
            @info "Post-processing DirectFlatPlate..."
            if backend == DOCKER
                work = "/tmp/DFP_post"
                foam_exec(backend, case_dir,
                    "rm -rf $work && cp -r . $work && cd $work" *
                    " && rm -rf dynamicCode postProcessing" *
                    " && simpleFoam -postProcess -time \"\$(ls -1d [0-9]* | grep -v '^0\$' | sort -g | tail -1)\"")
                pp_dest = joinpath(case_dir, "postProcessing")
                rm(pp_dest; force=true, recursive=true)
                run(`docker cp $(DOCKER_CONTAINER):$work/postProcessing $pp_dest`)
                foam_exec(backend, case_dir, "rm -rf $work")
            else
                foam_exec(backend, case_dir,
                    "rm -rf dynamicCode postProcessing" *
                    " && simpleFoam -postProcess -time \"\$(ls -1d [0-9]* | grep -v '^0\$' | sort -g | tail -1)\"")
            end
        end,

        :viz => () -> begin
            res    = plot_residuals(case_dir; savedir=plotting_dir, label="DirectFlatPlate")
            fields = plot_fields(case_dir; savedir=plotting_dir)
            wval   = nothing
            if inp.VAL.valPlot
                wval = plot_dfp_w_validation(case_dir; savedir=plotting_dir, gen=inp.VAL.Gen)
            end

            return (residuals=res, fields=fields, wval=wval)
        end,

        :_order => () -> [:clean, :prep, :mesh, :run, :post, :viz],
    )
end

MODULE_REGISTRY["DirectFlatPlate"] = make_direct_flat_plate
