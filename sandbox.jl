#!/usr/bin/env julia
# sandbox.jl - A testing ground for learning how SequenceVariation.jl works

# TODO: Try to variant call from BAM records
# TODO: Try to replicate `Haplotype` methods with `Variant` methods
# TODO: Try to export `Variation`s to VCF

using BioAlignments
using BioGenerics
using BioSequences
using BioSymbols
using FASTX
using GenomicFeatures: Strand
using SequenceVariation
using Statistics
using XAM

edits(v::Variant) = v.edits
Base.:(==)(x::Variation, y::Variation) = x.ref == y.ref && x.edit == y.edit
Base.hash(x::Variation, h::UInt) = hash(Variation, hash((x.ref, x.edit), h))

function variations(v::Variant)
    variations = Vector{Variation}(undef, length(v.edits))
    for (i, e) in enumerate(v.edits)
        variations[i] = Variation(v.ref, e)
    end #for
    return variations
end #function

function variations(vs::AbstractVector{Variant})
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

function BioGenerics.leftposition(v::Variation)
    return v.edit.pos
end #function

function readpos(v::Variation, a::Union{Alignment,AlignedSequence,PairwiseAlignment})
    return first(ref2seq(a, leftposition(v)))
end #function

function relativepos(v::Variation, b::BAM.Record)
    if leftposition(v) <= BAM.position(b)
        return 0.0
    elseif leftposition(v) >= BAM.rightposition(b)
        return 1.0
    else
        return readpos(v, BAM.alignment(b)) / (BAM.rightposition(b) - BAM.position(b))
    end #if
end #function

function substitution_quality(v::Variation, b::BAM.Record)
    return BAM.quality(b)[first(ref2seq(BAM.alignment(b), leftposition(v)))]
end #function

function insertion_quality(v::Variation, b::BAM.Record)
    if leftposition(v) <= BAM.position(b)
        startpos = first(ref2seq(BAM.alignment(b), BAM.position(b)))
        endpos = first(ref2seq(BAM.alignment(b), leftposition(v) + 1))
    elseif leftposition(v) >= BAM.rightposition(b)
        startpos = first(ref2seq(BAM.alignment(b), leftposition(v) - 1))
        endpos = first(ref2seq(BAM.alignment(b), BAM.rightposition(b)))
    else
        startpos = first(ref2seq(BAM.alignment(b), leftposition(v)))
        endpos = first(ref2seq(BAM.alignment(b), leftposition(v) + 1)) - 1
    end #if

    return mean(BAM.quality(b)[startpos:endpos])
end #function

function deletion_quality(v::Variation, b::BAM.Record)
    leftpos = first(ref2seq(BAM.alignment(b), leftposition(v) - 1))
    rightpos = first(ref2seq(BAM.alignment(b), leftposition(v) + length(v.edit.x)))

    return mean(BAM.quality(b)[leftpos:rightpos])
end #function

function quality(v::Variation, b::BAM.Record)
    if v.edit.x isa Substitution
        return substitution_quality(v, b)
    elseif v.edit.x isa Insertion
        return insertion_quality(v, b)
    elseif v.edit.x isa Deletion
        return deletion_quality(v, b)
    else
        return NaN
    end #if
end #function

struct VariationInfo{S<:BioSequence,T<:BioSymbol}
    variation::Variation{S,T}
    readpos::Float64
    quality::Float64
    strand::Strand
end #struct

variation(vi::VariationInfo) = vi.variation
readpos(vi::VariationInfo) = vi.readpos
quality(vi::VariationInfo) = vi.quality
strand(vi::VariationInfo) = vi.strand

function annotated_variations(query::BAM.Record, reference::NucleotideSeq)
    aligned_seq = AlignedSequence(BAM.sequence(query), BAM.alignment(query))
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
                BAM.ispositivestrand(query) ? Strand('+') : Strand('-'),
            ),
        )
    end #for
    return query_variations
end #function

read01 = AlignedSequence(
    dna"TTTATCTGTGTGAACTTCTTGGCTTAGTTT", Alignment("15M2P15M12H", 1, 6)
)
read02 = AlignedSequence(dna"CTGTGTGAACTTCTTGGCTTAGTATCGTTG", Alignment("30M", 1, 11))
refseq = dna"ACAACTTTATCTCTCTCAACTTCTTCCCTTACTATCCTTCACAACAATCCACACATTACTGCACTTTAAACACTTTTTTA"

aln01 = PairwiseAlignment(read01, refseq)
aln02 = PairwiseAlignment(read02, refseq)

var01 = Variant(aln01)
var02 = Variant(aln02)

bam_file = "example/sample.bam"
bai_file = "example/sample.bam.bai"
ref_file = "example/reference.fasta"

ref_records = FASTA.Record[]
FASTA.Reader(open(ref_file, "r")) do reader
    global ref_records = collect(reader)
end #do

reference_names = FASTA.identifier.(ref_records)

all_variants = VariationInfo[]

record = BAM.Record()
open(BAM.Reader, bam_file; index=bai_file) do reader
    while !eof(reader)
        empty!(record)
        read!(reader, record)

        if BAM.ismapped(record)
            if BAM.refname(record) in reference_names
                ref_strand = FASTA.sequence(
                    first(
                        filter(r -> FASTA.identifier(r) == BAM.refname(record), ref_records)
                    ),
                )
                push!(all_variants, annotated_variations(record, ref_strand)...)
            end #if
        end #if
    end #while
end #do

variation_pos_pileup = Dict{Variation,Vector{Float64}}()
variation_quality_pileup = Dict{Variation,Vector{Float64}}()
variation_strand_pileup = Dict{Variation,Vector{Strand}}()

variation_pos_avg = Dict{Variation,Float64}()
variation_quality_avg = Dict{Variation,Float64}()
variation_strand_avg = Dict{Variation,Float64}()
variation_altdepth = Dict{Variation,Int}()

for var in all_variants
    variation_pos_pileup[variation(var)] = Float64[]
    variation_quality_pileup[variation(var)] = Float64[]
    variation_strand_pileup[variation(var)] = Strand[]
end #for

for var in all_variants
    push!(variation_pos_pileup[variation(var)], readpos(var))
    push!(variation_quality_pileup[variation(var)], quality(var))
    push!(variation_strand_pileup[variation(var)], strand(var))
end #for

for (var, poss) in variation_pos_pileup
    variation_altdepth[var] = length(poss)
    variation_pos_avg[var] = mean(poss)
end #for

for (var, quals) in variation_quality_pileup
    variation_quality_avg[var] = mean(quals)
end #for

for (var, strd) in variation_strand_pileup
    n_pos = count(s -> s == Strand('+'), strd)
    N = length(strd)

    variation_strand_avg[var] = n_pos / N
end #for
