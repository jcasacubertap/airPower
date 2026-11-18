#
# airPower — load the central configuration into `inp`
#
# inputs.jl (`inputs`, set per case/run) and settings.jl (`settings`, set once,
# rarely touched) share one block layout and are merged block by block into the
# single `inp` every stage reads, so code accesses inp.<CaseName>.<variable>
# whichever file the entry lives in. A key must live in exactly one of the two:
# the same leaf in both is an error, never a silent override.
#
# Include once per process, from anywhere: the airPower root is found from this
# file's own location.
#

const AIRPOWER_ROOT_INPUTS = dirname(dirname(dirname(@__DIR__)))   # .../airPower

include(joinpath(AIRPOWER_ROOT_INPUTS, "inputs.jl"))     # const inputs
let f = joinpath(AIRPOWER_ROOT_INPUTS, "settings.jl")
    isfile(f) || error("settings.jl not found next to inputs.jl in $(AIRPOWER_ROOT_INPUTS). " *
                       "It holds the secondary entries split off inputs.jl; copy it from the " *
                       "canonical airPower.")
end
include(joinpath(AIRPOWER_ROOT_INPUTS, "settings.jl"))   # const settings

# Block-wise merge: keys of `a` in their order, then the keys only in `b`. A key
# in both is merged further when both sides are blocks (NamedTuples); a leaf in
# both is an error.
function mergeInputs(a::NamedTuple, b::NamedTuple, path::String = "inp")
    ks = vcat(collect(keys(a)), [k for k in keys(b) if !haskey(a, k)])
    vals = map(ks) do k
        if haskey(a, k) && haskey(b, k)
            va, vb = a[k], b[k]
            (va isa NamedTuple && vb isa NamedTuple) ||
                error("$(path).$(k) is set in both inputs.jl and settings.jl; keep it in one.")
            mergeInputs(va, vb, "$(path).$(k)")
        else
            haskey(a, k) ? a[k] : b[k]
        end
    end
    return NamedTuple{Tuple(ks)}(Tuple(vals))
end

const inp = mergeInputs(inputs, settings)
