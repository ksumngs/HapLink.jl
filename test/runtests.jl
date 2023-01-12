using BioAlignments:
    EDNAFULL,
    AlignedSequence,
    PairwiseAlignment,
    AffineGapScoreModel,
    GlobalAlignment,
    cigar,
    pairalign
using BioGenerics: leftposition
using BioSequences: BioSequence, LongDNA, @dna_str, ungap!
using FASTX: FASTA
using GenomicFeatures: Interval, Strand, STRAND_POS, STRAND_NEG
using Random: randstring, seed!
using SequenceVariation: Haplotype, Variation, variations
using Statistics: mean
using XAM: SAM, BAM
using VariantCallFormat: VCF

using Documenter: DocMeta, doctest
using Test
using TestItems
using TestItemRunner

using HapLink

seed!(1)

const DNA_MODEL = AffineGapScoreModel(EDNAFULL; gap_open=-12, gap_extend=-3)
align(a::BioSequence, b::BioSequence) = pairalign(GlobalAlignment(), a, b, DNA_MODEL).aln
function sam(a::PairwiseAlignment)
    qname = randstring()
    flag = 0x00
    rname = "REFERENCE"
    pos = a.a.aln.firstref
    mapq = 255
    cig = cigar(a.a.aln)
    rnext = "*"
    pnext = 0
    tlen = 0
    seq = a.a.seq
    qual = String(Char.(0x22:(0x22 + UInt8(length(seq)))))
    tag = "Xs:Z:HapLink"
    return SAM.Record(
        join([qname, flag, rname, pos, mapq, cig, rnext, pnext, tlen, seq, qual, tag], "\t")
    )
end #function
function write_sam(recs::AbstractVector{SAM.Record})
    sam_file, _ = mktemp()
    sam_header = SAM.Header([SAM.MetaInfo("SQ", ["SN" => "REFERENCE", "LN" => 35])])
    sam_writer = SAM.Writer(open(sam_file, "w"), sam_header)
    for rec in recs
        write(sam_writer, rec)
    end #for
    close(sam_writer)
    return sam_file
end #function

# Sequences from Japanese Encephalitis Virus full genome positions NC_001437:5979-7016 with
# an additional two-base deletion in genotype III
#! format: off
const REFERENCE    = ungap!(dna"CATCAGGGCTGACTGGA---TTGCCAAGCATGGCACT")
const GENOTYPE_I   = ungap!(dna"CATCAGGACTGACCGGA---TTGCCAAGCATGGCACT")
const GENOTYPE_II  = ungap!(dna"CACCAGGATTGACTGGA---TTGCCAA--ATGGCGCT")
const GENOTYPE_III = ungap!(dna"CATCAGGACTGACTGGA---TTGCCAAGCATGGCACT")
const GENOTYPE_V   = ungap!(dna"CATCCAGCGTGCCTGGAAGTCTGTCAAGCCTGGCGCT")
#! format: on

const GENOTYPES = [REFERENCE, GENOTYPE_I, GENOTYPE_II, GENOTYPE_III, GENOTYPE_V]
const ALIGNMENTS = align.(GENOTYPES, [REFERENCE])
const VARIANTS = Haplotype.(ALIGNMENTS)
const VARIATIONS = unique!(variations(VARIANTS))
const SAMS = sam.(ALIGNMENTS)
const SAM_FILE = write_sam(SAMS)
const VARIATIONINFOS = VariationInfo.(VARIATIONS, [0.5], [30.0], [STRAND_POS])
const VARIATIONPILEUPS =
    VariationPileup.(
        VARIATIONS, [0x04], [[0.0, 1.0]], [[20.0, 30.0]], [[STRAND_POS, STRAND_NEG]]
    )
const VARIATIONCALLS = call_variant.(VARIATIONPILEUPS, [1.0])
const VCFS = vcf.(VARIATIONCALLS, ["REFERENCE"])

@testset "HapLink.jl" begin
    include("doctests.jl")
    include("fasta.jl")
    include("xam.jl")
    include("variation.jl")
    include("variationinfo.jl")
    include("variationpileup.jl")
    include("variationcall.jl")
    include("findset.jl")
    include("haplotype.jl")
    @run_package_tests
end #testset
