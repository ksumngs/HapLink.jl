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

"""
    _posin(ps::Pseudoread, v::Variation)

Determines if the positions that make up `v` are wholly contained by `ps`
"""
function _posin(ps::Pseudoread, v::Variation)
    return leftposition(ps) <= leftposition(v) && rightposition(ps) >= rightposition(v)
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
