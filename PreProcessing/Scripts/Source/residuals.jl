using Glob

"""
    plot_residuals(case_path; savedir, label)

Parse `postProcessing/solverInfo/*/solverInfo.dat` and plot convergence history.
Saves PNG to `savedir/residuals_<label>.png`.
"""
function plot_residuals(case_path::AbstractString;
                        savedir::AbstractString=case_path,
                        label::AbstractString=basename(case_path))
    search_dir = joinpath(case_path, "postProcessing", "solverInfo")
    if !isdir(search_dir)
        @warn "No solverInfo directory in $case_path"
        return nothing
    end
    files = glob("*/solverInfo.dat", search_dir)
    if isempty(files)
        @warn "No solverInfo.dat found in $case_path"
        return nothing
    end

    # Use the first (or only) solverInfo.dat
    dat_file = first(files)
    @info "Parsing residuals from $dat_file"

    lines = readlines(dat_file)

    # Parse header (second comment line has column names)
    header_line = ""
    for l in lines
        if startswith(l, "# Time") || startswith(l, "#  Time")
            header_line = l
        end
    end
    if isempty(header_line)
        @warn "Could not find header in $dat_file"
        return nothing
    end

    cols = split(replace(header_line, "#" => ""))
    col_names = [strip(c) for c in cols if !isempty(strip(c))]

    # Find columns ending with "_initial" (these are the initial residuals)
    residual_cols = [(i, name) for (i, name) in enumerate(col_names) if endswith(name, "_initial")]

    # Parse data lines
    data_lines = [l for l in lines if !startswith(l, "#") && !isempty(strip(l))]
    if isempty(data_lines)
        @warn "No data lines in $dat_file"
        return nothing
    end

    nrows = length(data_lines)
    ncols = length(col_names)
    data = zeros(nrows, ncols)
    for (r, line) in enumerate(data_lines)
        fields = split(line)
        for (c, f) in enumerate(fields)
            if c <= ncols
                val = tryparse(Float64, f)
                data[r, c] = isnothing(val) ? NaN : val
            end
        end
    end

    iterations = data[:, 1]

    p = plot(; xlabel="Iteration", ylabel="Initial residual",
             title="Residuals — $label", yscale=:log10, legend=:topright,
             size=(800, 500), linewidth=1.5, grid=true)

    for (col_idx, col_name) in residual_cols
        field_name = replace(col_name, "_initial" => "")
        plot!(p, iterations, data[:, col_idx]; label=field_name)
    end

    mkpath(savedir)
    outfile = joinpath(savedir, "residuals_$(label).png")
    savefig(p, outfile)
    @info "Saved: $outfile"
    return p
end
