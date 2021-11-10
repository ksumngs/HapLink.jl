include("variant.jl")

export Haplotype

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
