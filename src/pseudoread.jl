struct Pseudoread
    startpos::UInt64
    endpos::UInt64
    read::Variant
end #struct

function Pseudoread(startpos::Integer, endpos::Integer, read::Variant)
    (startpos < 0 || endpos < 0) && error("Positions cannot be less than zero")
    return Pseudoread(UInt64(startpos), UInt64(endpos), read)
end #function

function Pseudoread(query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq)
    return Pseudoread(
        _XAM_(query).position(query),
        _XAM_(query).rightposition(query),
        Variant(query, reference),
    )
end #function

BioGenerics.leftposition(ps::Pseudoread) = ps.startpos
BioGenerics.rightposition(ps::Pseudoread) = ps.endpos
variant(ps::Pseudoread) = ps.read
_variations(ps::Pseudoread) = variations(variant(ps))

"""
    pseudoreads(sam::Union{AbstractString,AbstractPath}, consensus::NucleotideSeq)

Create a set of [`Pseudoread`](@ref)s from the alignments present in `sam` when aligned
against `consensus`.
"""
function pseudoreads(sam::Union{AbstractString,AbstractPath}, consensus::NucleotideSeq)
    XAM = _issam(sam) ? SAM : BAM

    returned_reads = Pseudoread[]

    reader = open(XAM.Reader, sam)
    for r in reader
        _is_primary_record(r) && push!(returned_reads, Pseudoread(r, consensus))
    end #for
    close(reader)

    return returned_reads
end #function

function simulate(
    pool::AbstractArray{Pseudoread},
    refseq::BioSequence,
    subcon::AbstractArray{Variation{S,T}};
    reverse_order::Union{Bool,Nothing}=nothing,
    overlap_min::Int64=0,
    overlap_max::Int64=500,
) where {S,T}
    # Get the reference sequence to start with
    refvar = Variant(refseq, Variation{S,T}[])

    # Find out if we're going to walk the genome forward or backward
    is_reversed = isnothing(reverse_order) ? rand(Bool) : reverse_order

    # Create a pseudoread that we will build from. It needs to start at the end of the
    # genome that we're starting the simulation from
    previous_read = if is_reversed
        Pseudoread(length(refseq), length(refseq), refvar)
    else
        Pseudoread(1, 1, refvar)
    end

    # Set up the order to go through variations
    var_pool = is_reversed ? reverse(sort(subcon)) : sort(subcon)

    first_pass = true

    # Add a read for every variation position
    for var in var_pool
        # Get reads that contain `leftpostion(var)` and `rightposition(var)`
        # and that overlap `previous_read` by the specified amount
        # and that match all existing `variations(previous_read)`
        # Overlaps do not need to be checked if this is the first read
        matching_read_pool = filter(
            r -> _is_candidate(
                r,
                previous_read,
                var,
                overlap_min,
                overlap_max;
                check_overlap=!first_pass,
            ),
            pool,
        )

        # Abort if we've hit a dead end
        if isempty(matching_read_pool)
            return missing
        end #if

        # Get the new read
        matching_read = rand(matching_read_pool)

        # Get the variations we've worked with so far as a vector
        previous_variations = _variations(previous_read)

        # Add any new variations to the old ones
        new_variations = union(_variations(matching_read), var_pool)
        push!(previous_variations, new_variations...)
        unique!(previous_variations)
        sort!(previous_variations)

        # Find out what the total interval we've transversed is
        start_pos = is_reversed ? leftposition(matching_read) : leftposition(previous_read)
        end_pos = is_reversed ? rightposition(previous_read) : rightposition(matching_read)

        # Create the new pseudoread
        previous_read = Pseudoread(start_pos, end_pos, Variant(refseq, previous_variations))

        first_pass = false
    end #for

    return variant(previous_read)
end #function

function _is_candidate(
    current_read::Pseudoread,
    previous_read::Pseudoread,
    var::Variation,
    overlap_min::Integer,
    overlap_max::Integer;
    check_overlap::Bool=true,
)
    _posin(current_read, var) || return false
    if check_overlap
        overlap_inrange(previous_read, current_read, overlap_min, overlap_max) ||
            return false
    end
    variations_match(previous_read, current_read) || return false
    return true
end #function

"""
    _posin(ps::Pseudoread, v::Variation)

Determines if the positions that make up `v` are wholly contained by `ps`
"""
function _posin(ps::Pseudoread, v::Variation)
    return leftposition(ps) <= leftposition(v) && rightposition(ps) >= rightposition(v)
end #function

"""
    variations_match(reference::Pseudoread, query::Union{Pseudoread,Haplotype})

Returns `true` if `query` could be a read from `reference`.

Additional `Variation`s found within `query` that are not present in `reference` **do not**
invalid a positive match result so long as none of those `Variation`s [`contridict`](@ref)
`reference`. These `Variation`s can either

1. Extend beyond the length of `reference`, or
2. Exist as "sequencing errors" within the interval of `reference`

On the other hand, **every** `Variation` within the overlapping interval of `reference` and
`query` that is present in `reference` **must** also be found in `query`. In other words,
`query` must have all of `reference` (within the overlap), but `reference` does not need all
of `query`.

This arrangment allows for the creation of bodies of matching [`Pseudoread`](@ref)s, as
done via [`simulate`](@ref).
"""
function variations_match(reference::Pseudoread, query::Haplotype)
    # Every variation within the reference must also be found within the query, provided
    # that the query covers that interval
    all(var -> var in query, overlapping_variations(reference, query)) || return false

    # Every variation in query can either...
    # 1. Exactly match one found in reference
    # 2. Not contridict one found in reference
    for var in variations(query)
        var in variant(reference) && continue

        contridicts(var, variant(reference)) && return false
    end #for

    return true
end #function

function variations_match(reference::Pseudoread, query::Pseudoread)
    return variations_match(reference, variant(query))
end #function

"""
    overlapping_variations(ps::Pseudoread, v::Variant)

Find all `Variation`s within `v` that are contained within the range defined by `ps`
"""
function overlapping_variations(ps::Pseudoread, v::Variant)
    return filter(
        var ->
            leftposition(var) >= leftposition(ps) &&
                rightposition(var) <= rightposition(ps),
        variations(v),
    )
end #function

"""
    overlap_inrange(
        x::Pseudoread, y::Pseudoread, overlap_min::Integer, overlap_max::Integer
    )

Returns `true` if the overlap between `x` and `y` is within `overlap_min` and `overlap_max`,
returns `false` otherwise

# Example
```jldoctest
julia> using BioSequences, SequenceVariation

julia> ps1 = Pseudoread(1, 100, Variant(dna"A", Variation{LongDNA{4}, DNA}[]));

julia> ps2 = Pseudoread(50, 150, Variant(dna"A", Variation{LongDNA{4}, DNA}[]));

julia> overlap_inrange(ps1, ps2, 1, 500)
true

julia> overlap_inrange(ps1, ps2, 75, 500)
false

julia> overlap_inrange(ps1, ps2, 1, 25)
false
```
"""
function overlap_inrange(
    x::Pseudoread, y::Pseudoread, overlap_min::Integer, overlap_max::Integer
)
    x_interval = leftposition(x):rightposition(x)
    y_interval = leftposition(y):rightposition(y)

    overlap_interval = overlap(x_interval, y_interval)
    overlap_amount = last(overlap_interval) - first(overlap_interval)

    return overlap_amount >= overlap_min && overlap_amount <= overlap_max
end #function

"""
    overlap(x::UnitRange{S}, y::UnitRange{S}) where {S}

Finds the inclusive overlap interval of `x` and `y`

# Example
```jldoctest
julia> overlap(1:5, 3:10)
3:5

julia> overlap(2:4, 6:8)
6:5

julia> overlap(1:10, 2:9)
2:9
```
"""
function overlap(x::UnitRange{S}, y::UnitRange{S}) where {S}
    a = first(x)
    b = first(y)
    c = last(x)
    d = last(y)

    e = max(a, b)
    f = min(c, d)

    return UnitRange{S}(e, f)
end #function
