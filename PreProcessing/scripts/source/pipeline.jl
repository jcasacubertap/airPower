const DEFAULT_ALL_ACTIONS = [:clean, :prep, :mesh, :solve, :post, :viz]

# Registry populated after module files are included
const MODULE_REGISTRY = Dict{String, Function}()

"""
    run_module(name, action_str; root)

Main entry point. `name` is the module name (e.g. "DirectFlatPlate"),
`action_str` is one of "clean", "mesh", "solve", "post", "viz", or "all".
`root` is the airPower project root directory.
"""
function run_module(name::AbstractString, action_str::AbstractString="all";
                    root::AbstractString)
    if !haskey(MODULE_REGISTRY, name)
        available = join(keys(MODULE_REGISTRY), ", ")
        error("Unknown module '$name'. Available: $available")
    end

    backend = detect_backend()
    @info "Backend: $backend"

    mod = MODULE_REGISTRY[name](backend, root)

    actions = if action_str == "all"
        haskey(mod, :_order) ? mod[:_order]() : DEFAULT_ALL_ACTIONS
    else
        [Symbol(action_str)]
    end

    result = nothing
    for action in actions
        if !haskey(mod, action)
            avail = filter(k -> k != :_order, collect(keys(mod)))
            error("Module '$name' does not support action :$action. " *
                  "Available: $(join(avail, ", "))")
        end
        @info "=== $name : $action ==="
        result = mod[action]()
    end
    @info "Done."
    return result
end
