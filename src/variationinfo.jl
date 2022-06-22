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

Calls `Variation`s based on the alignment in `query` against `reference`, and returns every
variation call found within `query` as a `Vector{VariationInfo}`
"""
@generated function variationinfos(
    query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq
)
    XAM = _xam_switch(query)

    quote
        aligned_seq = AlignedSequence($XAM.sequence(query), $XAM.alignment(query))
        paired_alignment = PairwiseAlignment(aligned_seq, reference)

        query_variant = Variant(paired_alignment)

        query_variations = VariationInfo[]

        for variation in variations(query_variant)
            push!(
                query_variations,
                VariationInfo(
                    variation,
                    relativepos(variation, query),
                    Float64(quality(variation, query)),
                    # Note: ispositivestrand is implemented for BAM records, but not for
                    # SAM records, so inline the logic here
                    $XAM.flag(query) & 0x10 == 0 ? Strand('+') : Strand('-'),
                ),
            )
        end #for
        return query_variations
    end #quote
end #function
