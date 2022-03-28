export Haplotype
export mutations

"""
    Haplotype(mutations::AbstractVector{Variant})

Create an object to describe a group of mutations that appear together. A single
[`Variant`](@ref) can be passed to create a haplotype of a single variant.
"""
struct Haplotype
    mutations::AbstractVector{Variant}
end #struct

function Haplotype(var::Variant)
    return Haplotype([var])
end #function

function Haplotype(hapdict::Dict{String,Any})
    v = Variant.(hapdict["mutations"])
    return Haplotype(v)
end

function Haplotype(::Nothing)
    return Haplotype(Variant[])
end #function

function Base.show(io::IO, h::Haplotype)
    return print(
        io,
        string(
            "Haplotype (",
            length(h.mutations),
            ") [",
            join([string(v.chromosome, ":", v.position) for v in h.mutations], ","),
            "]",
        ),
    )
end #function

function Base.isless(h1::Haplotype, h2::Haplotype)
    # Haplotypes with more variants are larger in magnitude
    # Otherwise rank based on the first variant position
    return length(h1.mutations) < length(h2.mutations) &&
           first(sort(h1.mutations)) < first(sort(h2.mutations))
end #function

function mutations(h::Haplotype)
    return h.mutations
end #function

function serialize_yaml(h::Haplotype; reason::Union{String,Nothing}=nothing)
    return string(
        "---\n",
        isnothing(reason) ? "" : string("reason: ", reason, "\n"),
        "mutations: \n",
        serialize_yaml.(h.mutations)...,
    )
end
