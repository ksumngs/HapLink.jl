function SequenceVariation.Haplotype(
    query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq
)
    XAM = typeof(query) == SAM.Record ? SAM : BAM

    aligned_seq = AlignedSequence(XAM.sequence(query), XAM.alignment(query))
    paired_alignment = PairwiseAlignment(aligned_seq, reference)

    return Haplotype(paired_alignment)
end #function

"""
    cigar(hap::Haplotype{S,T}) where {S,T}

Constructs a CIGAR string representing the alignment of the sequence of `hap` to its
reference.
"""
function cigar(hap::Haplotype{S,T}) where {S,T}
    cigar_string = String[]

    mismatch_vars = filter(var -> !isa(mutation(var), Substitution), variations(hap))

    length(mismatch_vars) > 0 || return "$(length(reference(hap)))M"

    lastvar = first(mismatch_vars)

    leftposition(lastvar) > 1 && push!(cigar_string, "$(leftposition(lastvar))M")

    for var in mismatch_vars
        push!(cigar_string, _cigar_between(lastvar, var))
        push!(cigar_string, _cigar(var))
        lastvar = var
    end #for

    remaining_bases = length(reference(hap)) - rightposition(lastvar)
    remaining_bases > 0 && push!(cigar_string, "$(remaining_bases)M")

    return join(cigar_string, "")
end #function
