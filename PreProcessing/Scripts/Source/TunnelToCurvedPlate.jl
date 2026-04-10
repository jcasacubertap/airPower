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
            if backend == DOCKER
                # VirtioFS workaround: write to /tmp, copy back
                foam_exec(backend, tunnel_case,
                          "surfaceOrient constant/triSurface/$(base).stl '$(outside_pt)' /tmp/$(base)_oriented.stl" *
                          " && cp /tmp/$(base)_oriented.stl constant/triSurface/$(base)_oriented.stl")
            else
                foam_exec(backend, tunnel_case,
                          "surfaceOrient constant/triSurface/$(base).stl '$(outside_pt)' constant/triSurface/$(base)_oriented.stl")
            end

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
            stl_cmd = "surfaceTransformPoints" *
                      " -scale '($chord $chord 1)'" *
                      " -translate '($dx 0 $dy)'" *
                      " -origin '($xo 0 $yo)'" *
                      " -rollPitchYaw '(0 0 $(-aoa))'" *
                      " constant/triSurface/$(base)_oriented.stl"
            if backend == DOCKER
                foam_exec(backend, tunnel_case, stl_cmd *
                          " /tmp/airfoil.stl && cp /tmp/airfoil.stl constant/triSurface/airfoil.stl")
            else
                foam_exec(backend, tunnel_case, stl_cmd * " constant/triSurface/airfoil.stl")
            end

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
            @info "Meshing TunnelCase..."
            if backend == DOCKER
                work = "/tmp/TunnelCase_mesh"
                foam_exec(backend, tunnel_case,
                    "rm -rf $work && cp -r . $work && cd $work" *
                    " && blockMesh && snappyHexMesh")
                tc_docker = docker_case_path(tunnel_case)
                for d in ["constant/polyMesh", "constant/triSurface", "constant/extendedFeatureEdgeMesh", "1", "2"]
                    run(ignorestatus(`docker exec $(DOCKER_CONTAINER) bash -c
                        "[ -d $work/$d ] && rm -rf $tc_docker/$d && cp -r $work/$d $tc_docker/$d"`))
                end
                foam_exec(backend, tunnel_case, "rm -rf $work")
            else
                foam_exec(backend, tunnel_case, "blockMesh")
                foam_exec(backend, tunnel_case, "snappyHexMesh")
            end
            @info "Meshing complete"
        end,

        :meshAirfoil => () -> begin
            write_airfoil_le_input_param(airfoil_case)

            @info "Generating AirfoilLECase grid..."
            run_julia_subprocess(generate_grid; dir=airfoil_case)
            if backend == DOCKER
                # Push blockMeshDict into Docker (VirtioFS truncates large files)
                bmd = joinpath(airfoil_case, "system", "blockMeshDict")
                run(`docker cp $bmd $(DOCKER_CONTAINER):$(docker_case_path(airfoil_case))/system/blockMeshDict`)
            end
            @info "Meshing AirfoilLECase..."
            foam_exec(backend, airfoil_case, "blockMesh")
        end,

        :runTunnel => () -> begin
            write_tunnel_input_param(tunnel_case)
            @info "Solving TunnelCase..."
            if backend == DOCKER
                work = "/tmp/TunnelCase_run"
                np = inp.TTCP.nProcs
                foam_exec(backend, tunnel_case,
                    "rm -rf $work && cp -r . $work && cd $work" *
                    " && cp -r 2 2.bak" *
                    " && cp 0/U 0/p 0/k 0/omega 0/nut 2/" *
                    " && decomposePar -force -time 2" *
                    " && mpirun -np $np --oversubscribe simpleFoam -parallel" *
                    " && reconstructPar" *
                    " && rm -r 2 && mv 2.bak 2")
                latest = read(`docker exec $(DOCKER_CONTAINER) bash -c
                    "ls -1d $work/[0-9]* 2>/dev/null | grep -v '^$work/0\$' | grep -v '^$work/2\$' | sort -g | tail -1"`, String) |> strip
                if !isempty(latest)
                    run(`docker cp $(DOCKER_CONTAINER):$latest $(tunnel_case)/`)
                    @info "Copied time directory: $(basename(latest))"
                end
                foam_exec(backend, tunnel_case, "rm -rf $work")
            else
                foam_script(backend, tunnel_case, "run", "$(inp.TTCP.nProcs)")
            end
        end,

        :map => () -> begin
            @info "Mapping tunnel → airfoil..."
            run_julia_subprocess(map_script; dir=tunnel_to_curved)
        end,

        :runAirfoil => () -> begin
            write_airfoil_le_input_param(airfoil_case)
            @info "Solving AirfoilLECase..."
            if backend == DOCKER
                work = "/tmp/AirfoilLECase_run"
                foam_exec(backend, airfoil_case,
                    "rm -rf $work && cp -r . $work && cd $work" *
                    " && decomposePar -force" *
                    " && mpirun -np $(inp.TTCP.nProcs) --oversubscribe simpleFoam -parallel" *
                    " && reconstructPar")
                latest = read(`docker exec $(DOCKER_CONTAINER) bash -c
                    "ls -1d $work/[0-9]* 2>/dev/null | sort -g | tail -1"`, String) |> strip
                if !isempty(latest) && basename(latest) != "0"
                    run(`docker cp $(DOCKER_CONTAINER):$latest $(airfoil_case)/`)
                    @info "Copied time directory: $(basename(latest))"
                end
                foam_exec(backend, airfoil_case, "rm -rf $work")
            else
                foam_script(backend, airfoil_case, "run", "$(inp.TTCP.nProcs)")
            end
        end,

        :postAirfoil => () -> begin
            @info "Post-processing AirfoilLECase..."
            if backend == DOCKER
                work = "/tmp/AirfoilLECase_post"
                foam_exec(backend, airfoil_case,
                    "rm -rf $work && cp -r . $work && cd $work" *
                    " && rm -rf dynamicCode postProcessing" *
                    " && simpleFoam -postProcess -time \"\$(ls -1d [0-9]* | sort -g | tail -1)\"")
                pp_dest = joinpath(airfoil_case, "postProcessing")
                rm(pp_dest; force=true, recursive=true)
                run(`docker cp $(DOCKER_CONTAINER):$work/postProcessing $pp_dest`)
                foam_exec(backend, airfoil_case, "rm -rf $work")
            else
                foam_exec(backend, airfoil_case,
                    "rm -rf dynamicCode postProcessing" *
                    " && simpleFoam -postProcess -time \"\$(ls -1d [0-9]* | sort -g | tail -1)\"")
            end
            @info "Post-processing complete"
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
            a = inp.TTCP.airfoilLE
            geom_args = (chord_mm=t.chord * 1000, alpha_deg=t.alphaDeg,
                         x_center_mm=t.xCenter * 1000, y_center_mm=t.yCenter * 1000)
            fields = plot_fields(airfoil_case; savedir=plotting_dir,
                delta=0.010, geom_args...)
            wall = plot_wall_geometry(airfoil_case; savedir=plotting_dir, geom_args...)
            wallq = plot_wall_quantities(airfoil_case; savedir=plotting_dir, geom_args...)
            expval = nothing
            if inp.VAL.valPlot
                expval = plot_experimental_validation(airfoil_case;
                    savedir=plotting_dir, gen=inp.VAL.Gen, delta=0.010, geom_args...)
            end

            return (res_airfoil=res_airfoil, fields=fields, wall=wall, wallq=wallq,
                    expval=expval)
        end,

        :monitorTunnel => () -> begin
            if backend == DOCKER
                tmp_dir = mktempdir()
                dst = joinpath(tmp_dir, "postProcessing", "solverInfo", "0")
                mkpath(dst)
                for f in ["solverInfo_0.dat", "solverInfo.dat"]
                    src = "/tmp/TunnelCase_run/postProcessing/solverInfo/0/$f"
                    run(ignorestatus(`docker cp $(DOCKER_CONTAINER):$src $dst/`))
                end
                f0 = joinpath(dst, "solverInfo_0.dat")
                fd = joinpath(dst, "solverInfo.dat")
                if isfile(f0) && (!isfile(fd) || filesize(fd) < filesize(f0))
                    cp(f0, fd; force=true)
                end
                p = plot_residuals(tmp_dir; savedir=tmp_dir, label="TunnelCase (live)")
                rm(tmp_dir; force=true, recursive=true)
            else
                p = plot_residuals(tunnel_case; savedir=mktempdir(), label="TunnelCase (live)")
            end
            display(p)
            return p
        end,

        :monitorAirfoil => () -> begin
            if backend == DOCKER
                tmp_dir = mktempdir()
                dst = joinpath(tmp_dir, "postProcessing", "solverInfo", "0")
                mkpath(dst)
                for f in ["solverInfo_0.dat", "solverInfo.dat"]
                    src = "/tmp/AirfoilLECase_run/postProcessing/solverInfo/0/$f"
                    run(ignorestatus(`docker cp $(DOCKER_CONTAINER):$src $dst/`))
                end
                f0 = joinpath(dst, "solverInfo_0.dat")
                fd = joinpath(dst, "solverInfo.dat")
                if isfile(f0) && (!isfile(fd) || filesize(fd) < filesize(f0))
                    cp(f0, fd; force=true)
                end
                p = plot_residuals(tmp_dir; savedir=tmp_dir, label="AirfoilLECase (live)")
                rm(tmp_dir; force=true, recursive=true)
            else
                p = plot_residuals(airfoil_case; savedir=mktempdir(), label="AirfoilLECase (live)")
            end
            display(p)
            return p
        end,

        :_order => () -> [:clean, :prep, :meshTunnel, :runTunnel,
                          :meshAirfoil, :map, :runAirfoil, :postAirfoil,
                          :vizTunnel, :vizAirfoil],
    )
end

MODULE_REGISTRY["TunnelToCurvedPlate"] = make_tunnel_to_curved_plate
