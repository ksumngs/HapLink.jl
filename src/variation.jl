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
    XAM = _xam_record_switch(r)

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

@generated function _subqual(v::Variation, r::Union{SAM.Record,BAM.Record})
    XAM = _xam_record_switch(r)

    quote
        # Substitution quality: basecall quality of substituted base
        # ref    GATTACA
        #        ||| |||
        # seq    GATAACA  => qscore = '%' = 4
        # qscore "#$%&'(
        return $XAM.quality(r)[first(ref2seq($XAM.alignment(r), leftposition(v)))]
    end #quote
end #function

@generated function _insqual(v::Variation, r::Union{SAM.Record,BAM.Record})
    XAM = _xam_record_switch(r)

    quote
        if leftposition(v) <= $XAM.position(r)
            # Insertion quality if insertion starts before alignment: basecall quality of
            # read from first base until last inserted base
            # ref    --GATTACA
            #          |||||||
            # seq    GGGATTACA => qscore = mean('#$') = mean([1,2]) = 1.5
            # qscore "#$%&'()*
            startpos = 1
            endpos = first(ref2seq($XAM.alignment(r), leftposition(v))) - 1
        elseif leftposition(v) >= $XAM.rightposition(r)
            startpos = first(ref2seq($XAM.alignment(r), leftposition(v) - 1)) + 1
            endpos = length($XAM.quality(r))
        else
            # Insertion quality: basecall quality of inserted bases
            # ref    GAT--TACA
            #        |||  ||||
            # seq    GATGGTACA  => qscore = mean('%&') = mean([4,5]) = 4.5
            # qscore "#$%&'()*
            startpos = first(ref2seq($XAM.alignment(r), leftposition(v) - 1)) + 1
            endpos = first(ref2seq($XAM.alignment(r), leftposition(v))) - 1
        end #if

        return mean($XAM.quality(r)[startpos:endpos])
    end #quote
end #function

@generated function _delqual(v::Variation, r::Union{SAM.Record,BAM.Record})
    XAM = _xam_record_switch(r)

    quote
        leftpos = first(ref2seq($XAM.alignment(r), leftposition(v) - 1))
        rightpos = first(ref2seq($XAM.alignment(r), leftposition(v) + length(mutation(v))))

        return mean($XAM.quality(r)[leftpos:rightpos])
    end #quote
end #function

"""
    quality(v::Variation, r::Union{SAM.Record,BAM.Record}) -> Float64

Get the phred-scalled basecall quality of `v` within the sequencing read of `r`.
"""
function quality(v::Variation, r::Union{SAM.Record,BAM.Record})
    if mutation(v) isa Substitution
        return _subqual(v, r)
    elseif mutation(v) isa Insertion
        return _insqual(v, r)
    elseif mutation(v) isa Deletion
        return _delqual(v, r)
    else
        return NaN
    end #if
end #function

function interval(v::Variation, refname::AbstractString)
    return Interval(refname, Int(leftposition(v)), Int(rightposition(v)))
end #function
