"""
    SequenceVariation.Haplotype(
        query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq
    ) -> SequenceVariation.Haplotype

Specialized constructor that allows converting the alignment in a `XAM.Record` into a
`Haplotype`.
"""
function SequenceVariation.Haplotype(
    query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq
)
    XAM = typeof(query) == SAM.Record ? SAM : BAM

    aligned_seq = AlignedSequence(XAM.sequence(query), XAM.alignment(query))
    paired_alignment = PairwiseAlignment(aligned_seq, reference)

    return Haplotype(paired_alignment)
end #function

function SequenceVariation.Haplotype(ref::NucleotideSeq, data::AbstractDict)
    snps = data["snps"]
    vars = map(snp -> Variation(ref, snp), snps)
    return Haplotype(ref, vars)
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

function _name(h::Haplotype; prefix::AbstractString="", is_consensus::Bool=false)
    refseq = copy(reference(h))
    mutseq = reconstruct(h)

    seqhash = is_consensus ? "CONSENSUS" : _short_hash(mutseq)

    if isempty(prefix)
        return seqhash
    else
        return join([prefix, seqhash], "_")
    end #if
end #function

function _short_hash(x::Any; length::Integer=8)
    hash = bytes2hex(sha1(string(x)))
    return hash[1:length]
end #function
