#!/usr/bin/env julia
# sandbox.jl - A testing ground for learning how SequenceVariation.jl works

# TODO: Try to variant call from BAM records
# TODO: Try to replicate `Haplotype` methods with `Variant` methods
# TODO: Try to export `Variation`s to VCF

using BioAlignments
using FASTX
using SequenceVariation
using XAM

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

println(all_variants)
