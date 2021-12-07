using BioAlignments
using BioSequences
using BioSymbols
using XAM

export baseatreferenceposition
export matchvariant
export mutate

"""
    myref2seq(aln::Alignment, i::Int)

Replicates the functionality of BioAlignments `ref2seq`, but can handle hard clips
by effectively removing them for the intent of finding the position.
"""
function myref2seq(aln::Alignment, i::Int)
    if aln.anchors[2].op == OP_HARD_CLIP
        # Hard clipping was shown on operation 2
        # (operation 1 is always a start position)

        # Save where the clipping ends
        alnstart = aln.anchors[2]

        # Declare a new empty array where we can rebuild the alignment
        newanchors = AlignmentAnchor[]

        # Rebase the start of our new alignment to where the clipping ends
        push!(newanchors, AlignmentAnchor(0, aln.anchors[1].refpos, OP_START))

        # Add new anchors
        for j in 3:(length(aln.anchors) - 1)
            newanchor = AlignmentAnchor(
                aln.anchors[j].seqpos - alnstart.seqpos,
                aln.anchors[j].refpos,
                aln.anchors[j].op,
            )
            push!(newanchors, newanchor)
        end #for

        # Package up our new alignment
        workingalignment = Alignment(newanchors)
    else
        # Package up the old alignment if there was no hard clipping
        workingalignment = aln
    end #if

    # Check that the requested base is in range
    if !seqisinrange(workingalignment, i)
        return (0, OP_HARD_CLIP)
    end

    # Perform regular alignment search, minus any hard clipping
    return ref2seq(workingalignment, i)
end #function

function seqisinrange(aln::Alignment, i::Int)
    reflen = i - first(aln.anchors).refpos
    seqlen = last(aln.anchors).seqpos - first(aln.anchors).seqpos
    return seqlen > reflen
end #function

function firstseqpos(aln::Alignment)
    return first(aln.anchors).seqpos
end #function

function lastseqpos(aln::Alignment)
    return last(aln.anchors).seqpos
end #function

"""
    matchvariant(base::Union{NucleotideSeq,DNA,AbstractVector{DNA}}, var::Variant)

Checks if `base` matches the reference or variant expected in `var`, and returns a symbol
indicating which, if any, it matches.

Returned values can be `:reference` for a reference match, `:alternate` for an alternate
match, or `:other` for no match with the given variant.
"""
function matchvariant(base::NucleotideSeq, var::Variant)
    refbase = LongDNASeq(var.referencebase)
    altbase = LongDNASeq(var.alternatebase)

    if base == refbase
        return :reference
    elseif base == altbase
        return :alternate
    else
        return :other
    end #if
end #function

function matchvariant(base::DNA, var::Variant)
    return matchvariant(LongDNASeq([base]), var)
end

function matchvariant(base::AbstractVector{DNA}, var::Variant)
    return matchvariant(LongDNASeq(base), var)
end

"""
    baseatreferenceposition(record::BAM.Record, pos::Int)

Get the base at reference position `pos` present in the sequence of `record`.
"""
function baseatreferenceposition(record::BAM.Record, pos::Int)
    seqpos = myref2seq(BAM.alignment(record), pos)[1]
    if seqpos > 0 && seqpos < BAM.seqlength(record)
        return BAM.sequence(record)[seqpos]
    else
        return DNA_N
    end
end # function

"""
    function mutate(seq::FASTA.Record, haplotype::Haplotype)
    function mutate(seq::NucleotideSeq, haplotype::Haplotype)

Give the reference sequence `seq` mutations to match the position and basecalls of
`haplotype`. Returns a new sequence, leaving `seq` unmodified.

When mutating a `FASTA.Record`, the new record is given a new unique identifier and
description based on the SHA1 hash of the complete genotype.
"""
function mutate(record::FASTA.Record, haplotype::Haplotype)
    newseq = mutate(sequence(record), haplotype)
    sequencehash = bytes2hex(sha1(string(newseq)))
    newid = sequencehash[1:8]
    newdesc = string(description(record), ", variant ", sequencehash)

    return FASTA.Record(newid, newdesc, newseq)
end #function

function mutate(seq::NucleotideSeq, haplotype::Haplotype)
    newseq = seq

    for var in haplotype.mutations
        newseq[var.position] = var.alternatebase[1]
    end #for

    return newseq
end #function
