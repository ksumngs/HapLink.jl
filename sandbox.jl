#!/usr/bin/env julia
# sandbox.jl - A testing ground for learning how SequenceVariation.jl works

# TODO: Try to variant call from BAM records
# TODO: Try to replicate `Haplotype` methods with `Variant` methods
# TODO: Try to export `Variation`s to VCF

using BioAlignments
using BioSequences
using BioSymbols
using FASTX
using SequenceVariation
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

read01 = AlignedSequence(dna"TTTATCTGTGTGAACTTCTTGGCTTAGTTT", Alignment("30M", 1, 6))
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

all_variants = Variant[]

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

                aligned_seq = AlignedSequence(BAM.sequence(record), BAM.alignment(record))
                paired_alignment = PairwiseAlignment(aligned_seq, ref_strand)

                push!(all_variants, Variant(paired_alignment))
            end #if
        end #if
    end #while
end #do

all_variations = cat(variations.(all_variants)...; dims=1)
unique_variations = unique(all_variations)

for v in unique_variations
    containing_variants = filter(var -> v in var, all_variants)
    if !isempty(containing_variants)
        if v.edit.x isa Substitution{<:BioSymbol}
            pos = v.edit.pos
            ref = v.ref
            altdepth = count(var -> v in var, all_variants)
            println(
                "reference\t$pos\t.\t$(ref[pos])\t$(v.edit.x.x)\t.\tPASS\tAD=$(altdepth)"
            )
        end #if
    end #if
end #for
