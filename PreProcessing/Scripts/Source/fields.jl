using DelimitedFiles

"""
    plot_fields(case_path; savedir, filename)

Parse `postProcessing/midPlane.csv` and plot 2D contours of velocity magnitude and pressure.
CSV columns: x, y, z, u, v, w, p
"""
function plot_fields(case_path::AbstractString;
                     savedir::AbstractString=case_path,
                     filename::AbstractString="midPlane.csv")
    csv_path = joinpath(case_path, "postProcessing", filename)
    if !isfile(csv_path)
        @warn "midPlane.csv not found at $csv_path"
        return nothing
    end

    @info "Parsing fields from $csv_path"

    # Read CSV line-by-line, skipping header and malformed lines
    lines = readlines(csv_path)
    x = Float64[]; y = Float64[]; u = Float64[]; v = Float64[]
    w = Float64[]; p = Float64[]
    for line in lines[2:end]    # skip header
        fields = split(line, ',')
        length(fields) == 7 || continue
        vals = tryparse.(Float64, fields)
        any(isnothing, vals) && continue
        push!(x, vals[1]); push!(y, vals[2])
        push!(u, vals[4]); push!(v, vals[5]); push!(w, vals[6])
        push!(p, vals[7])
    end
    if isempty(x)
        @warn "No valid data rows in $csv_path"
        return nothing
    end
    @info "Parsed $(length(x)) cells"

    Umag = @. sqrt(u^2 + v^2 + w^2)

    mkpath(savedir)
    label = basename(case_path)

    # Auto-detect aspect ratio: use :equal for compact domains, :auto for thin ones
    dx = maximum(x) - minimum(x)
    dy = maximum(y) - minimum(y)
    ar = (dx > 0 && dy > 0 && dx / dy > 5) ? :auto : :equal

    # Velocity magnitude contour
    p1 = scatter(x, y; marker_z=Umag, markersize=1, markerstrokewidth=0,
                 aspect_ratio=ar, colorbar=true, colorbar_title="|U| [m/s]",
                 xlabel="x [m]", ylabel="y [m]",
                 title="Velocity magnitude — $label",
                 size=(1000, 600), legend=false)
    outfile1 = joinpath(savedir, "Umag_$(label).png")
    savefig(p1, outfile1)
    @info "Saved: $outfile1"

    # Pressure contour
    p2 = scatter(x, y; marker_z=p, markersize=1, markerstrokewidth=0,
                 aspect_ratio=ar, colorbar=true, colorbar_title="p [Pa]",
                 xlabel="x [m]", ylabel="y [m]",
                 title="Pressure — $label",
                 size=(1000, 600), legend=false)
    outfile2 = joinpath(savedir, "pressure_$(label).png")
    savefig(p2, outfile2)
    @info "Saved: $outfile2"

    return (p1, p2)
end
