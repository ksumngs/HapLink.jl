using BioAlignments:
    EDNAFULL,
    AlignedSequence,
    PairwiseAlignment,
    AffineGapScoreModel,
    GlobalAlignment,
    cigar,
    pairalign
using BioGenerics: leftposition
using BioSequences: BioSequence, @rna_str, ungap!
using Documenter: DocMeta, doctest
using HapLink
using SequenceVariation: Variant, Variation, variations
using Test
using XAM: SAM, BAM

const DNA_MODEL = AffineGapScoreModel(EDNAFULL; gap_open=-12, gap_extend=-3)
align(a::BioSequence, b::BioSequence) = pairalign(GlobalAlignment(), a, b, DNA_MODEL).aln
function sam(a::PairwiseAlignment)
    return SAM.Record(
        "*\t0\tREFERENCE\t$(a.a.aln.firstref)\t255\t$(cigar(a.a.aln))\t*\t0\t0\t$(a.a.seq)\t*\tXs:Z:HapLink",
    )
end #function

# Sequences from Japanese Encephalitis Virus full genome positions NC_001437:5979-7016 with
# an additional two-base deletion in genotype III
#! format: off
const REFERENCE    = ungap!(rna"CAUCAGGGCUGACUGGA---UUGCCAAGCAUGGCACU")
const GENOTYPE_I   = ungap!(rna"CAUCAGGACUGACCGGA---UUGCCAAGCAUGGCACU")
const GENOTYPE_II  = ungap!(rna"CACCAGGAUUGACUGGA---UUGCCAA--AUGGCGCU")
const GENOTYPE_III = ungap!(rna"CAUCAGGACUGACUGGA---UUGCCAAGCAUGGCACU")
const GENOTYPE_V   = ungap!(rna"CAUCCAGCGUGCCUGGAAGUCUGUCAAGCCUGGCGCU")
#! format: on

const GENOTYPES = [REFERENCE, GENOTYPE_I, GENOTYPE_II, GENOTYPE_III, GENOTYPE_V]
const ALIGNMENTS = align.(GENOTYPES, [REFERENCE])
const VARIANTS = Variant.(ALIGNMENTS)
const VARIATIONS = unique!(variations(VARIANTS))
const SAMS = sam.(ALIGNMENTS)

@testset "HapLink.jl" begin
    include("doctests.jl")
    include("variation.jl")
end #testset
