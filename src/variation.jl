"""
    variations(vs::AbstractVector{Variant})

Extracts all `SequenceVariation.Variation`s from `vs`.
"""
function SequenceVariation.variations(vs::AbstractVector{<:Variant})
    nvar = sum([length(v.edits) for v in vs])
    all_variations = Vector{Variation}(undef, nvar)
    i = 1

    for v in vs
        for var in variations(v)
            all_variations[i] = var
            i += 1
        end #for
    end #for
    return all_variations
end #function

"""
    seqpos(v::Variation, a::Union{Alignment,AlignedSequence,PairwiseAlignment})

Get the position of `v` in the sequence of `a`.

# Example
```jldoctest
using BioAlignments, BioSequences, SequenceVariation
v = Variation(dna"AAAAA", "A3T")
a = Alignment("2=1X2=", 1, 1)
seqpos(v, a)

# output

3
```
"""
function seqpos(v::Variation, a::Union{Alignment,AlignedSequence,PairwiseAlignment})
    if typeof(mutation(v)) <: Union{Substitution,Deletion}
        # Substitutions and deletions consume reference, so their positions match the edit's
        # position
        return first(ref2seq(a, leftposition(v)))
    elseif mutation(v) isa Insertion
        # Insertions do not consume reference, so back up one base from the reference
        # position of the base, find the sequence position for that reference, then
        # add one to advance the sequence into the insertion
        return first(ref2seq(a, Int(leftposition(v)) - 1)) + 1
    else
        error("Unknown Variation edit type $(typeof(mutation(v)))")
    end #if
end #function

"""
    relativepos(v::Variation, r::Union{SAM.Record,BAM.Record})

Calculates the fractional position of `v` within the sequence of `r`. If `v` is
out-of-bounds of `r`, then will return `0` for positions before `r` and `1` for positions
after `r`.
"""
@generated function relativepos(v::Variation, r::Union{SAM.Record,BAM.Record})
    XAM = _xam_switch(r)

    quote
        if leftposition(v) <= $XAM.position(r)
            return 0.0
        elseif leftposition(v) >= $XAM.rightposition(r)
            return 1.0
        else
            pos = seqpos(v, $XAM.alignment(r)) / ($XAM.rightposition(r) - $XAM.position(r))
            return pos <= 1.0 ? pos : 1.0
        end #if
    end #quote
end #function
