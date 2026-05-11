"""
    _plot_wall_modulation_check(case_dir, wm, savedir)

Plot expected bump shape (solid line) vs actual mesh wall coordinates (symbols),
zoomed to the refinement region.
"""
function _plot_wall_modulation_check(case_dir::AbstractString, wm, savedir::AbstractString)
    # Compute bump extent
    if wm.shape == :sigmoidal
        bump_xs = wm.xStart
        bump_xe = wm.xEnd
    elseif wm.shape == :esn
        bump_xs, bump_xe = _esn_geometry(wm)
    end
    bumpL = bump_xe - bump_xs

    # Plot range: refinement box1 extent (bump ± 50% bump width)
    x_lo = max(0.0, bump_xs - bumpL / 2)
    x_hi = bump_xe + bumpL / 2

    # Expected shape (dense sampling)
    x_expected = range(x_lo, x_hi, length=2000)
    h_expected = [wall_bump(xi, wm) * 1000.0 for xi in x_expected]  # mm

    # Actual mesh wall: parse plate patch face centres
    x_mesh = Float64[]
    y_mesh = Float64[]
    try
        xm, ym = parse_patch_face_centers(case_dir, "plate")
        # Filter to plot range and one z-layer
        for i in eachindex(xm)
            if x_lo <= xm[i] <= x_hi
                push!(x_mesh, xm[i])
                push!(y_mesh, ym[i] * 1000.0)  # mm
            end
        end
    catch e
        @warn "Could not parse mesh wall: $e"
    end

    mkpath(savedir)

    common_opts = (
        xlabel     = "x [m]",
        ylabel     = "h [mm]",
        framestyle = :box,
        grid       = true,
        gridalpha  = 0.3,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 9,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        size       = (900, 400),
        dpi        = 200,
        legend     = :topright,
    )

    # Symbols first (underneath)
    if !isempty(x_mesh)
        fig = scatter(x_mesh, y_mesh;
            label      = "Mesh wall (face centres)",
            color      = :firebrick,
            markersize = 2,
            markerstrokewidth = 0,
            common_opts...)
    else
        fig = plot(; common_opts...)
    end

    # Solid line on top
    plot!(fig, x_expected, h_expected;
        label      = "Expected shape",
        color      = :black,
        linewidth  = 2,
    )

    outfile = joinpath(savedir, "wallModulationCheck.png")
    savefig(fig, outfile)
    @info "Saved: $outfile"
end

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

            wm = merge(inp.wallModulation, inp.DFP.wallBump)
            refine_cmds = ""
            # Note: topoSet/refineMesh disabled — it distorts the bump geometry
            # (face centres are averaged, not re-evaluated on the polyLine).
            # Use gridXfactor for resolution instead.

            if backend == DOCKER
                work = "/tmp/DFP_mesh"
                foam_exec(backend, case_dir,
                    "rm -rf $work && cp -r . $work && cd $work" *
                    " && blockMesh" * refine_cmds)
                dfp_docker = docker_case_path(case_dir)
                run(ignorestatus(`docker exec $(DOCKER_CONTAINER) bash -c
                    "rm -rf $dfp_docker/constant/polyMesh && cp -r $work/constant/polyMesh $dfp_docker/constant/polyMesh"`))
                foam_exec(backend, case_dir, "rm -rf $work")
            else
                foam_exec(backend, case_dir, "blockMesh" * refine_cmds)
            end

            # Plot bump vs actual mesh wall in the refinement region
            if wm.enabled
                @info "Generating wall modulation verification plot..."
                _plot_wall_modulation_check(case_dir, wm, plotting_dir)
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
            fields = plot_fields(case_dir; savedir=plotting_dir,
                                 wm=merge(inp.wallModulation, inp.DFP.wallBump))
            wval   = nothing
            if inp.VAL.valPlot
                wval = plot_dfp_w_validation(case_dir; savedir=plotting_dir, gen=inp.VAL.Gen, case_id=inp.VAL.Case)
            end

            return (residuals=res, fields=fields, wval=wval)
        end,

        :_order => () -> [:clean, :prep, :mesh, :run, :post, :viz],
    )
end

MODULE_REGISTRY["DirectFlatPlate"] = make_direct_flat_plate
