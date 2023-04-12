struct HaplotypeCall{S<:BioSequence,T<:BioSymbol}
    depth::UInt64
    linkage::Float64
    frequency::Float64
    significance::Float64
    haplotype::Haplotype{S,T}
end #struct

depth(hc::HaplotypeCall) = hc.depth
linkage(hc::HaplotypeCall) = hc.linkage
frequency(hc::HaplotypeCall) = hc.frequency
significance(hc::HaplotypeCall) = hc.significance
haplotype(hc::HaplotypeCall) = hc.haplotype
variant(hc::HaplotypeCall) = haplotype(hc)

function HaplotypeCall(
    haplotype::Haplotype{S,T}, reads::AbstractArray{Haplotype{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    hapcounts = occurence_matrix(haplotype, reads)

    depth = last(hapcounts)
    frequency = last(hapcounts) / sum(hapcounts)
    Δ = linkage(hapcounts)
    p = significance(hapcounts)

    return HaplotypeCall(depth, Δ, frequency, p, haplotype)
end #function

function HaplotypeCall(
    haplotype::AbstractArray{Variation{S,T}}, reads::AbstractArray{Haplotype{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    refseq = reference(first(haplotype))
    hapvar = Haplotype(refseq, haplotype)

    return HaplotypeCall(hapvar, reads)
end #function

function _name(hc::HaplotypeCall; prefix::AbstractString="")
    refseq = copy(reference(variant(hc)))
    mutseq = reconstruct(variant(hc))

    seqhash = bytes2hex(sha1(string(mutseq)))
    shorthash = seqhash[1:8]

    if isempty(prefix)
        return shorthash
    else
        return join([prefix, shorthash], "_")
    end #if
end #function

"""
    function ishaplotype(
        haplotype::Union{AbstractArray{Variation{S,T}}, Haplotype{S,T}},
        reads::AbstractArray{Haplotype{S,T}};
        frequency::Union{Float64,Nothing}=nothing,
        significance::Union{Float64,Nothing}=nothing,
        depth::Union{Int,Nothing}=nothing,
    ) where {S<:BioSequence,T<:BioSymbol}

Determines if a call of `haplotype` is supported by the sequences in `reads` based upon
the provided keyword criteria.

# Arguments
- `haplotype::Union{AbstractArray{Variation}, Haplotype}`: A `Vector` of `Variation`s or a
    `Haplotype` to search for as a haplotype
- `reads::AbstractArray{Haplotype}`: The reads to search for `haplotype` in

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
    reads::AbstractArray{Haplotype{S,T}};
    min_frequency::Union{Float64,Nothing}=nothing,
    significance_level::Union{Float64,Nothing}=nothing,
    min_depth::Union{Int,Nothing}=nothing,
) where {S<:BioSequence,T<:BioSymbol}
    check_contridicts(haplotype) && return false
    hap = Haplotype(reference(first(haplotype)), haplotype)

    return ishaplotype(
        hap,
        reads;
        min_frequency=min_frequency,
        significance_level=significance_level,
        min_depth=min_depth,
    )
end #function

function ishaplotype(
    haplotype::Haplotype{S,T},
    reads::AbstractArray{Haplotype{S,T}};
    min_frequency::Union{Float64,Nothing}=nothing,
    significance_level::Union{Float64,Nothing}=nothing,
    min_depth::Union{Int,Nothing}=nothing,
) where {S<:BioSequence,T<:BioSymbol}
    call = HaplotypeCall(haplotype, reads)

    if !isnothing(min_depth)
        depth(call) >= min_depth || return false
    end #if

    if !isnothing(min_frequency)
        frequency(call) >= min_frequency || return false
    end #if

    if !isnothing(significance_level)
        significance(call) <= significance_level || return false
    end #if

    return true
end #function

function check_contridicts(vars::AbstractArray{Variation{S,T}}) where {S,T}
    passed_vars = similar(vars, 0)
    ref = reference(first(vars))

    for var in vars
        passed_hap = Haplotype(ref, passed_vars)
        contridicts(var, passed_hap) && return true
        push!(passed_vars, var)
    end #for

    return false
end #function

function _dict(hc::HaplotypeCall; prefix::AbstractString="")
    variants = String[]
    for var in variations(variant(hc))
        varbuffer = IOBuffer()
        show(varbuffer, var)
        push!(variants, String(take!(varbuffer)))
        close(varbuffer)
    end #for
    return OrderedDict{String,Any}(
        "name" => _name(hc; prefix=prefix),
        "depth" => depth(hc),
        "frequency" => frequency(hc),
        "linkage" => linkage(hc),
        "significance" => significance(hc),
        "snps" => variants,
    )
end #function

"""
    occurence_matrix(
        haplotype::AbstractArray{Variation{S,T}},
        reads::AbstractArray{Haplotype{S,T}},
    ) where {S<:BioSequence,T<:BioSymbol}
    occurence_matrix(
        haplotype::Haplotype{S,T},
        reads::AbstractArray{S,T}
    ) where {S<:BioSequence,T<:BioSymbol}

Determine how many times the variants in `haplotype` appear in `reads` as an ``N``-
dimensional matrix.

# Arguments
- `haplotype::AbstractArray{Variation}`: A `Vector` of `Variation`s or a `Haplotype` to search
    for as a haplotype
- `reads::AbstractArray{Haplotype}`: The reads to search for `haplotype` in

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
    haplotype::AbstractArray{Variation{S,T}}, reads::AbstractArray{Haplotype{S,T}}
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
    haplotype::Haplotype{S,T}, reads::AbstractArray{Haplotype{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    return occurence_matrix(variations(haplotype), reads)
end #function

"""
    linkage(counts::AbstractArray{<:Integer})

Calculates the linkage disequilibrium of a haplotype given its ``N``-dimensional contingency
matrix, `counts`.

`counts` is an ``N``-dimensional array where the ``N``th dimension represents the ``N``th
variant call locus within a haplotype. `findoccurrences` produces such an array.

# Extended help

`linkage(::AbstractArray{<:Integer})` calculates an unweighted linkage disequilibrium as
given by Equation (6) of [Slatkin (1972)](https://doi.org/10.1093/genetics/72.1.157).

```math
D(1..N) = E\\left( \\prod_k^N i_k - P_k \\right)
```

where

- ``N`` is the number of variant loci
- ``D(1..N)`` is the linkage disequilibrium of all ``N`` variant loci
- ``E`` is an operator returning the arithmetic mean of its argument over every read
- ``i_k`` is a function that returns ``1`` if the ``k``-th locus of the given read contains
  the reference allele, and returns ``0`` otherwise.
- ``P_k`` is the probability of any read containing the reference allele in the ``k``-th
  locus, i.e. the frequency at which the reference allele is found within the entire read set
  at the ``k``-th locus
"""
function linkage(counts::AbstractArray{<:Integer})
    # Make math easier by declaring variables
    N = ndims(counts)
    n = sum(counts)

    # Short-circuit if there is only one contingency
    N > 1 || return 0.0

    return Σ(map(i -> _e_operand(i, counts) * counts[i], CartesianIndices(counts))) / n
end #function

function _e_operand(coord::CartesianIndex, counts::AbstractArray{<:Integer})
    # Skip if there are no reads for this contingency
    counts[coord] > 0 || return 0

    # Make math nicer-looking
    N = ndims(counts)

    # Return the result of the operand
    return Π(j -> _sub_i(j, coord, counts), 1:N)
end #function

function _e_square_operand(loci::Integer, counts::AbstractArray{<:Integer})
    return sum(i -> _sub_i(loci, i, counts)^2 * counts[i], CartesianIndices(counts)) /
           sum(counts)
end #function

function _sub_i(
    loci_index::Integer, read_index::CartesianIndex, counts::AbstractArray{<:Integer}
)
    is_reference = read_index[loci_index] == 1

    return Int(is_reference) - P_ref(counts, loci_index)
end #function

"""
    significance(counts::AbstractArray{<:Integer})

Calculates the ``χ^2`` significance of a haplotype given its ``N``-dimensional contingency
matrix, `counts`

See the documentation on [`linkage(::AbstractArray{<:Integer})`](@ref) or
[`occurence_matrix`](@ref) for details on how `counts` should be constructed.
"""
function significance(counts::AbstractArray{<:Integer})
    Χ_squared = _correlation_coeff(counts)^2 * sum(counts)

    # Calculate the significance
    return 1 - cdf(Chisq(1), Χ_squared)
end #function

function _correlation_coeff(counts::AbstractArray{<:Integer})
    return linkage(counts) / √Π(j -> _e_square_operand(j, counts), 1:ndims(counts))
end #function

function P_ref(counts::AbstractArray, dim::Integer)
    return sumsliced(counts, dim) / sum(counts)
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
