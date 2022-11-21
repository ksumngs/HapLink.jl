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
function relativepos(v::Variation, r::Union{SAM.Record,BAM.Record})
    if leftposition(v) <= _XAM_(r).position(r)
        return 0.0
    elseif leftposition(v) >= _XAM_(r).rightposition(r)
        return 1.0
    else
        pos =
            seqpos(v, _XAM_(r).alignment(r)) /
            (_XAM_(r).rightposition(r) - _XAM_(r).position(r))
        return pos <= 1.0 ? pos : 1.0
    end #if
end #function

function _subqual(v::Variation, r::Union{SAM.Record,BAM.Record})
    # Substitution quality: basecall quality of substituted base
    # ref    GATTACA
    #        ||| |||
    # seq    GATAACA  => qscore = '%' = 4
    # qscore "#$%&'(
    return _XAM_(r).quality(r)[first(ref2seq(_XAM_(r).alignment(r), leftposition(v)))]
end #function

function _insqual(v::Variation, r::Union{SAM.Record,BAM.Record})
    if leftposition(v) <= _XAM_(r).position(r)
        # Insertion quality if insertion starts before alignment: basecall quality of
        # read from first base until last inserted base
        # ref    --GATTACA
        #          |||||||
        # seq    GGGATTACA => qscore = mean('#$') = mean([1,2]) = 1.5
        # qscore "#$%&'()*
        startpos = 1
        endpos = first(ref2seq(_XAM_(r).alignment(r), leftposition(v))) - 1
    elseif leftposition(v) >= _XAM_(r).rightposition(r)
        startpos = first(ref2seq(_XAM_(r).alignment(r), leftposition(v) - 1)) + 1
        endpos = length(_XAM_(r).quality(r))
    else
        # Insertion quality: basecall quality of inserted bases
        # ref    GAT--TACA
        #        |||  ||||
        # seq    GATGGTACA  => qscore = mean('%&') = mean([4,5]) = 4.5
        # qscore "#$%&'()*
        startpos = first(ref2seq(_XAM_(r).alignment(r), leftposition(v) - 1)) + 1
        endpos = first(ref2seq(_XAM_(r).alignment(r), leftposition(v))) - 1
    end #if

    return mean(_XAM_(r).quality(r)[startpos:endpos])
end #function

function _delqual(v::Variation, r::Union{SAM.Record,BAM.Record})
    leftpos = first(ref2seq(_XAM_(r).alignment(r), leftposition(v) - 1))
    rightpos = first(ref2seq(_XAM_(r).alignment(r), leftposition(v) + length(mutation(v))))

    return mean(_XAM_(r).quality(r)[leftpos:rightpos])
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

"""
    variation(r::VCF.Record, refseq::NucleotideSeq)

Construct a `Variation` from `r` applying to `refseq`. There is no validation that `r`'s
actually describes a mutation in `refseq`.
"""
function variation(r::VCF.Record, refseq::NucleotideSeq)
    # Convert VCF record to array of nucleotides
    SEQTYPE = typeof(refseq)
    altbases = SEQTYPE(first(VCF.alt(r)))
    refbases = SEQTYPE(VCF.ref(r))
    pos = VCF.pos(r)

    if length(altbases) == 1 && length(refbases) == 1
        # Substitution
        return Variation(refseq, "$refbases$pos$altbases")
    elseif length(altbases) > length(refbases)
        # Insertion
        totalalt = pos > 1 ? altbases[2:end] : altbases[1:(end - 1)]
        return Variation(refseq, "$pos$totalalt")
    elseif length(refbases) > length(altbases)
        dellen = length(refbases) - length(altbases)
        delendpos = pos + dellen - 1
        return Variation(refseq, "Δ$pos-$delendpos")
    end #if

    return altseq
end #function

"""
    subconsensus_variations(vcf::Union{AbstractPath,AbstractString}, consensus::Variant)

Get a `Vector{Variation}` with passing variant calls from `vcf` that do not appear in
`consensus`
"""
function subconsensus_variations(
    vcf::Union{AbstractPath,AbstractString}, consensus::Variant{S,T}
) where {S,T}
    subcon_vars = Variation{S,T}[]

    reference_sequence = reference(consensus)

    reader = VCF.Reader(open(string(vcf), "r"))
    for record in reader
        if all(f -> f == "PASS", VCF.filter(record))
            var = variation(record, reference_sequence)

            if !(var in consensus)
                push!(subcon_vars, var)
            end #if
        end #if
    end #for
    close(reader)

    return subcon_vars
end #function
