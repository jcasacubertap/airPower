#!/usr/bin/env julia
#
# seedInputs — copy the numerical input values of an ORIGIN case onto a
# DESTINATION case's inputs.jl and settings.jl, preserving the destination's
# structure, comments and any NEW entries it has introduced since the origin
# was split off.
#
# Usage:
#   julia PreProcessing/src/selfRunning/seedInputs.jl <originCase> <destinationCase>
#
#   <originCase>       case whose values you want to reuse (the variant).
#   <destinationCase>  case to update in place (the fresh canonical).
#
# Each argument may be a case directory (its inputs.jl, plus settings.jl when
# present) or a single .jl file. The configuration of a case is the union of its
# files, so an entry is matched by its full path (e.g. TTCP.tunnel.chord)
# whichever file holds it: an older single-file inputs.jl seeds both inputs.jl
# and settings.jl of a split destination. Every leaf present on both sides takes
# the origin's value verbatim. Entries only in the destination (new in canonical)
# keep their current value and are reported for review; entries only in the
# origin (since removed) are reported too. Each destination file is backed up to
# <file>.bak before writing.

# --- resolve a case dir or a .jl path to its configuration files ----------
# A case directory gives inputs.jl, plus settings.jl when it exists (cases made
# before the split have inputs.jl alone).
function resolve_inputs(path::AbstractString)
    if isdir(path)
        f = joinpath(path, "inputs.jl")
        isfile(f) || error("no inputs.jl in case directory: $(path)")
        s = joinpath(path, "settings.jl")
        return isfile(s) ? [f, s] : [f]
    elseif isfile(path)
        return [path]
    else
        error("not found: $(path)")
    end
end

# --- line classification --------------------------------------------------

# Strip a trailing line comment (# ...), respecting double-quoted strings.
function strip_comment(line::AbstractString)
    instr = false
    for (i, c) in pairs(line)
        if c == '"'
            instr = !instr
        elseif c == '#' && !instr
            return line[1:prevind(line, i)]
        end
    end
    return line
end

# Leaf assignments on a (comment-stripped) code fragment -> [(key, valueText)].
# Value is a quoted string, a :symbol, or a run up to the next top-level comma.
const ASSIGN_RE = r"([A-Za-z_]\w*)\s*=\s*(\"[^\"]*\"|:[A-Za-z_]\w*|[^,]+?)\s*(?:,|$)"

leaf_assignments(code::AbstractString) =
    [(String(m.captures[1]), strip(m.captures[2])) for m in eachmatch(ASSIGN_RE, code)]

# Classify a raw line -> (:root | :open | :close | :leaf | :skip, payload).
function classify(line::AbstractString)
    code = strip(strip_comment(line))
    isempty(code) && return (:skip, nothing)
    if occursin(r"^const\s+[A-Za-z_]\w*\s*=\s*\($", code)   # inp | inputs | settings
        return (:root, nothing)
    elseif (m = match(r"^([A-Za-z_]\w*)\s*=\s*\($", code)) !== nothing
        return (:open, String(m.captures[1]))
    elseif startswith(code, ")")
        return (:close, nothing)
    else
        return (:leaf, code)
    end
end

# --- value substitution on an original line (keeps key, spacing, comment) ---
function set_value(line::AbstractString, key::AbstractString, newval::AbstractString)
    pat = Regex("((?<![A-Za-z0-9_])" * key * "(?![A-Za-z0-9_])\\s*=\\s*)" *
                "(\"[^\"]*\"|:[A-Za-z_]\\w*|[^,#]+?)(?=\\s*(?:,|#|\$))")
    ok = false
    newline = replace(line, pat => function (s)
        ok = true
        String(match(pat, s).captures[1]) * newval
    end; count = 1)
    return (newline, ok)
end

# --- collect path -> valueText from a file's lines -------------------------
function collect_values(lines)
    stack = String[]
    vals  = Dict{Vector{String},String}()
    for line in lines
        (kind, payload) = classify(line)
        if kind === :root
            push!(stack, "__root__")
        elseif kind === :open
            push!(stack, payload)
        elseif kind === :close
            isempty(stack) || pop!(stack)
        elseif kind === :leaf
            base = filter(!=("__root__"), stack)
            for (k, v) in leaf_assignments(payload)
                vals[vcat(base, k)] = v
            end
        end
    end
    return vals
end

# --- rewrite destination lines, pulling matched values from the origin map --
function rewrite(lines, originMap)
    stack    = String[]
    out      = String[]
    copied   = 0
    destOnly = Vector{String}[]
    failed   = Vector{String}[]
    used     = Set{Vector{String}}()
    for line in lines
        (kind, payload) = classify(line)
        if kind === :root
            push!(stack, "__root__"); push!(out, line)
        elseif kind === :open
            push!(stack, payload); push!(out, line)
        elseif kind === :close
            isempty(stack) || pop!(stack); push!(out, line)
        elseif kind === :leaf
            base    = filter(!=("__root__"), stack)
            newline = line
            for (k, _) in leaf_assignments(payload)
                path = vcat(base, k)
                if haskey(originMap, path)
                    (newline, ok) = set_value(newline, k, originMap[path])
                    ok ? (copied += 1; push!(used, path)) : push!(failed, path)
                else
                    push!(destOnly, path)
                end
            end
            push!(out, newline)
        else
            push!(out, line)
        end
    end
    originOnly = [p for p in keys(originMap) if !(p in used)]
    return out, copied, destOnly, originOnly, failed
end

function main()
    length(ARGS) >= 2 ||
        error("usage: julia PreProcessing/src/selfRunning/seedInputs.jl <originCase> <destinationCase>")
    originFiles = resolve_inputs(ARGS[1])
    destFiles   = resolve_inputs(ARGS[2])
    isempty(intersect(abspath.(originFiles), abspath.(destFiles))) ||
        error("origin and destination share a file: nothing to seed")

    # the origin's configuration is the union of its files
    originMap = Dict{Vector{String},String}()
    for f in originFiles
        merge!(originMap, collect_values(readlines(f)))
    end

    copied   = 0
    destOnly = Vector{String}[]
    failed   = Vector{String}[]
    seeded   = Set{Vector{String}}()
    for f in destFiles
        out, n, dOnly, oOnly, fl = rewrite(readlines(f), originMap)
        copied += n
        append!(destOnly, dOnly);  append!(failed, fl)
        union!(seeded, setdiff(keys(originMap), oOnly))   # used in this file
        cp(f, f * ".bak"; force = true)
        open(f, "w") do io
            for l in out
                println(io, l)
            end
        end
    end
    originOnly = [p for p in keys(originMap) if !(p in seeded)]

    println("Seeded $(copied) value(s)")
    println("  origin:      ", join(originFiles, ", "))
    for f in destFiles
        println("  destination: $(f)   (backup: $(f).bak)")
    end
    if !isempty(destOnly)
        println("\n⚠  NEW in destination — kept current value, please review:")
        foreach(p -> println("     ", join(p, ".")), sort(destOnly))
    end
    if !isempty(originOnly)
        println("\nℹ  In origin only (obsolete / ignored):")
        foreach(p -> println("     ", join(p, ".")), sort(originOnly))
    end
    if !isempty(failed)
        println("\n✗  Matched but value not substituted (check by hand):")
        foreach(p -> println("     ", join(p, ".")), sort(failed))
    end
end

main()
