struct HaplotypeCall{S<:BioSequence,T<:BioSymbol}
    depth::UInt64
    linkage::Float64
    frequency::Float64
    significance::Float64
    haplotype::Variant{S,T}
end #struct

depth(hc::HaplotypeCall) = hc.depth
linkage(hc::HaplotypeCall) = hc.linkage
frequency(hc::HaplotypeCall) = hc.frequency
significance(hc::HaplotypeCall) = hc.significance
haplotype(hc::HaplotypeCall) = hc.haplotype
variant(hc::HaplotypeCall) = haplotype(hc)

function HaplotypeCall(
    haplotype::Variant{S,T}, reads::AbstractArray{Variant{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    hapcounts = occurence_matrix(haplotype, reads)

    depth = last(hapcounts)
    frequency = last(hapcounts) / sum(hapcounts)
    (Δ, p) = linkage(hapcounts)

    return HaplotypeCall(depth, Δ, frequency, p, haplotype)
end #function

function HaplotypeCall(
    haplotype::AbstractArray{Variation{S,T}}, reads::AbstractArray{Variant{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    refseq = reference(first(haplotype))
    hapvar = Variant(refseq, haplotype)

    return HaplotypeCall(hapvar, reads)
end #function

"""
    function ishaplotype(
        haplotype::Union{AbstractArray{Variation{S,T}}, Variant{S,T}},
        reads::AbstractArray{Variant{S,T}};
        frequency::Union{Float64,Nothing}=nothing,
        significance::Union{Float64,Nothing}=nothing,
        depth::Union{Int,Nothing}=nothing,
    ) where {S<:BioSequence,T<:BioSymbol}

Determines if a call of `haplotype` is supported by the sequences in `reads` based upon
the provided keyword criteria.

# Arguments
- `haplotype::Union{AbstractArray{Variation}, Variant}`: A `Vector` of `Variation`s or a
    `Variant` to search for as a haplotype
- `reads::AbstractArray{Variant}`: The reads to search for `haplotype` in

# Keywords
- `frequency::Union{Float64,Nothing}=nothing`: The minimum number of times the entire
    `haplotype` must appear within `reads` compared to the number of reads to return `true`
- `significance::Union{Float64,Nothing}=nothing`: The ``χ^2`` significance level (``α``)
    of linkage disequilibrium that a haplotype must surpass to return `true`
- `depth::Union{Int,Nothing}=nothing`: The minimum number of times the entire `haplotype`
    must appear within `reads` to return `true`

# Extended help
Linkage disequilibrium (``Δ``) is calculated by

```math
Δ = P_{reference} - \\prod_i P_{ref,i}
```

where
- ``P_{reference}`` is the probability of finding a read that contains only reference
    (consensus) bases
- ``P_{ref,i}`` is the probability of finding a read that contains the reference (consensus)
    base for variant ``i`` within a haplotype

and the ``χ^2`` statistic is calculated by

```math
χ^2 = \\frac{
            Δ^2 n
        }
        {
            \\left(\\prod_i^N {P_{ref,i} \\left(1 - P_{ref,i}\\right)}\\right)^\\frac{2}{N}
        }
```

where
- ``N`` is the number of `Variation`s in `haplotype` (i.e., `length(haplotype)`)
- ``n`` is the number of total reads sampled (i.e. `length(reads)`)

The significance is then calculated from the cumulative ``χ^2`` distribution function.
"""
function ishaplotype(
    haplotype::AbstractArray{Variation{S,T}},
    reads::AbstractArray{Variant{S,T}};
    frequency::Union{Float64,Nothing}=nothing,
    significance::Union{Float64,Nothing}=nothing,
    depth::Union{Int,Nothing}=nothing,
) where {S<:BioSequence,T<:BioSymbol}
    hapcounts = occurence_matrix(haplotype, reads)

    if !isnothing(depth)
        last(hapcounts) >= depth || return false
    end #if

    if !isnothing(frequency)
        last(hapcounts) / sum(hapcounts) >= frequency || return false
    end #if

    if !isnothing(significance)
        (Δ, p) = linkage(hapcounts)
        p <= significance || return false
    end #if

    return true
end #function

function ishaplotype(
    haplotype::Variant{S,T},
    reads::AbstractArray{Variant{S,T}};
    frequency::Union{Float64,Nothing}=nothing,
    significance::Union{Float64,Nothing}=nothing,
    depth::Union{Int,Nothing}=nothing,
) where {S<:BioSequence,T<:BioSymbol}
    return ishaplotype(
        variations(haplotype),
        reads;
        frequency=frequency,
        significance=significance,
        depth=depth,
    )
end #function

"""
    occurence_matrix(
        haplotype::AbstractArray{Variation{S,T}},
        reads::AbstractArray{Variant{S,T}},
    ) where {S<:BioSequence,T<:BioSymbol}
    occurence_matrix(
        haplotype::Variant{S,T},
        reads::AbstractArray{S,T}
    ) where {S<:BioSequence,T<:BioSymbol}

Determine how many times the variants in `haplotype` appear in `reads` as an ``N``-
dimensional matrix.

# Arguments
- `haplotype::AbstractArray{Variation}`: A `Vector` of `Variation`s or a `Variant` to search
    for as a haplotype
- `reads::AbstractArray{Variant}`: The reads to search for `haplotype` in

# Returns
- `2x2x... Array{Int, length(haplotype)}`: An ``N``-dimensional matrix where ``N`` is
    the number of variant positions in `readmatches`. The ``1`` index position in the
    ``n``th dimension represents the number of times the ``n``th variant position was found
    to have the reference base called, while the ``2`` index position represents the number
    of times the ``n``th variant position was found to have the alternate base called. E.g.
    `first(occurrence_matrix(reads))` gives the number of times the all-reference base
    haplotype was found in `reads`, while `occurrence_matrix(reads)[end]` gives the number
    of times the all-alternate base haplotype was found.
"""
function occurence_matrix(
    haplotype::AbstractArray{Variation{S,T}}, reads::AbstractArray{Variant{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    hapcounts = zeros(UInt, repeat([2], length(haplotype))...)

    for read in reads
        coordinates = zeros(Int, size(haplotype))
        i = 1

        for variation in haplotype
            if variation in read
                coordinates[i] = 1
            end #if
            i += 1
        end #for

        coordinates = CartesianIndex((coordinates .+ 1)...)
        hapcounts[coordinates] += 1
    end #for

    return hapcounts
end #function

function occurence_matrix(
    haplotype::Variant{S,T}, reads::AbstractArray{Variant{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    return occurence_matrix(variations(haplotype), reads)
end #function

"""
    linkage(counts::AbstractArray{Int})

Calculates the linkage disequilibrium and Chi-squared significance level of a combination of
haplotypes whose number of occurrences are given by `counts`.

`counts` is an ``N``-dimensional array where the ``N``th dimension represents the ``N``th
variant call position within a haplotype. `findoccurrences` produces such an array.
"""
function linkage(counts::AbstractArray{<:Integer})
    # Get the probability of finding a perfect reference sequence
    P_allref = first(counts) / sum(counts)

    # Get the probabilities of finding reference bases in any of the haplotypes
    P_refs = sumsliced.([counts], 1:ndims(counts)) ./ sum(counts)

    # Calculate linkage disequilibrium
    Δ = P_allref - prod(P_refs)

    # Calculate the test statistic
    r = Δ / (prod(P_refs .* (1 .- P_refs))^(1 / ndims(counts)))
    Χ_squared = r^2 * sum(counts)

    # Calculate the significance
    p = 1 - cdf(Chisq(1), Χ_squared)

    return Δ, p
end #function

"""
    sumsliced(A::AbstractArray, dim::Int, pos::Int=1)

Sum all elements that are that can be referenced by `pos` in the `dim` dimension of `A`.

# Example

```jldoctest
julia> A = reshape(1:8, 2, 2, 2)
2×2×2 reshape(::UnitRange{Int64}, 2, 2, 2) with eltype Int64:
[:, :, 1] =
 1  3
 2  4

[:, :, 2] =
 5  7
 6  8

julia> sumsliced(A, 2)
14

julia> sumsliced(A, 2, 2)
22
```

Heavily inspired by Holy, Tim "Multidimensional algorithms and iteration"
<https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes>
"""
function sumsliced(A::AbstractArray, dim::Int, pos::Int=1)
    i_pre = CartesianIndices(size(A)[1:(dim - 1)])
    i_post = CartesianIndices(size(A)[(dim + 1):end])
    return sum(A[i_pre, pos, i_post])
end #function
