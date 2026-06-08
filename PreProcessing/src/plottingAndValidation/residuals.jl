using Glob, LaTeXStrings

"""
    plot_residuals(case_path; savedir, label)

Parse `postProcessing/solverInfo/*/solverInfo.dat` and plot convergence history.
Saves PNG to `savedir/residuals<label>.png`.
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

    dat_file = first(files)
    @info "Parsing residuals from $dat_file"

    lines = readlines(dat_file)

    # Parse header
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

    # Map field names to LaTeX legend labels
    latex_names = Dict(
        "Ux" => L"U_x", "Uy" => L"U_y", "Uz" => L"U_z",
        "p"  => L"p",    "k"  => L"k",   "omega" => L"\omega",
    )

    colors = [:royalblue, :firebrick, :forestgreen, :darkorange, :purple, :teal]

    p = plot(;
        xlabel  = L"\mathrm{Iteration}",
        ylabel  = L"\mathrm{Initial \ residual}",
        yscale  = :log10,
        legend  = :topright,
        size    = (720, 480),
        framestyle   = :box,
        grid         = true,
        gridalpha    = 0.3,
        gridlinewidth = 0.5,
        minorgrid    = true,
        minorgridalpha = 0.15,
        tickfontsize   = 10,
        guidefontsize  = 12,
        legendfontsize = 9,
        titlefontsize  = 13,
        left_margin    = 8Plots.mm,
        bottom_margin  = 6Plots.mm,
        right_margin   = 4Plots.mm,
        top_margin     = 4Plots.mm,
        dpi = 200,
    )

    for (k, (col_idx, col_name)) in enumerate(residual_cols)
        field_name = replace(col_name, "_initial" => "")
        lbl = get(latex_names, field_name, latexstring(field_name))
        c = colors[mod1(k, length(colors))]
        plot!(p, iterations, data[:, col_idx];
              label=lbl, color=c, linewidth=2)
    end

    mkpath(savedir)
    outfile = joinpath(savedir, "residuals$(label).png")
    savefig(p, outfile)
    @info "Saved: $outfile"
    return p
end
