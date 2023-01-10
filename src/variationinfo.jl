"""
    VariationInfo{S<:BioSequence,T<:BioSymbol}

Represents statistics associated with a [`SequenceVariation.Variation`](@ref) within an
aligned sequencing read.

# Fields

- `variation::Variation{S,T}`: The `Variation` this object represents
- `readpos::Float64`: The location where `variation` occurs within a sequencing read
- `quality::Float64`: The phred-scaled basecall quality of `variation`
- `Strand::GenomicFeatures.Strand`: Which strand `variation` appears on
"""
struct VariationInfo{S<:BioSequence,T<:BioSymbol}
    variation::Variation{S,T}
    readpos::Float64
    quality::Float64
    strand::Strand
end #struct

"""
    variation(vi::VariationInfo) -> SequenceVariation.Variation

Gets the `Variation` of a [`VariationInfo`](@ref)
"""
variation(vi::VariationInfo) = vi.variation

"""
    readpos(vi::VariationInfo) -> Float64

Gets the position of `variation(vi)` within a sequencing read
"""
readpos(vi::VariationInfo) = vi.readpos

"""
    quality(vi::VariationInfo) -> Float64

Gets the phred-scaled basecall quality of `variation(vi)` within a sequencing read
"""
quality(vi::VariationInfo) = vi.quality

"""
    strand(vi::VariationInfo) -> GenomicFeatures.Strand

Gets the strand that `variation(vi)` appears on
"""
strand(vi::VariationInfo) = vi.strand

"""
    variationinfos(
        query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq
    ) -> Vector{VariationInfo}
    variationinfos(
        query::Union{AbstractString,AbstractPath},
        reference::Union{AbstractString,AbstractPath}
    ) -> Vector{VariationInfo}

Calls `Variation`s based on the alignments in `query` against `reference`, and returns every
variation call found within `query` as a `Vector{VariationInfo}`
"""
function variationinfos(query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq)
    query_variant = Haplotype(query, reference)

    query_variations = VariationInfo[]

    try
        for var in variations(query_variant)
            push!(
                query_variations,
                VariationInfo(
                    var,
                    relativepos(var, query),
                    Float64(quality(var, query)),
                    # Note: ispositivestrand is implemented for BAM records, but not for
                    # SAM records, so inline the logic here
                    _XAM_(query).flag(query) & 0x10 == 0 ? Strand('+') : Strand('-'),
                ),
            )
        end #for
    catch
        tname = _XAM_(query).tempname(query)
        @warn "Could not parse variations from $tname"
    end #try
    return query_variations
end #function

function variationinfos(
    query::Union{AbstractString,AbstractPath}, reference::Union{AbstractString,AbstractPath}
)
    refrecord = _first_record(reference)
    readtype = _issam(query) ? SAM.Reader : BAM.Reader
    reader = open(readtype, query)
    vi = _varinfos(reader, refrecord)
    return vi
end #struct

function _varinfos(reader::Union{SAM.Reader,BAM.Reader}, ref::FASTA.Record)
    all_variation_infos = VariationInfo[]
    record = _XAM_(reader).Record()

    while !eof(reader)
        empty!(record)
        read!(reader, record)

        if _is_primary_record(record) &&
            _XAM_(reader).refname(record) == FASTA.identifier(ref)
            seqtype = typeof(_XAM_(reader).sequence(record))
            push!(
                all_variation_infos, variationinfos(record, FASTA.sequence(seqtype, ref))...
            )
        end #if
    end #while

    return all_variation_infos
end #function

function _is_primary_record(r::Union{SAM.Record,BAM.Record})
    return _XAM_(r).ismapped(r) && (_XAM_(r).flag(r) & 0x900 == 0) # is primary line
end #function
