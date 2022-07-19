module HapLink

using ArgParse: ArgParseSettings, @add_arg_table!, project_version
using BioAlignments: Alignment, AlignedSequence, PairwiseAlignment, ref2seq
using BioGenerics: BioGenerics, leftposition, rightposition, metadata
using BioSequences: BioSequence, NucleotideSeq
using BioSymbols: BioSymbol
using Dates: Dates, today
using FASTX: FASTA
using FilePaths: FilePaths, AbstractPath, Path, absolute
using GenomicFeatures: Interval, Strand, STRAND_POS, eachoverlap
using HypothesisTests: FisherExactTest, pvalue
using SequenceVariation:
    SequenceVariation,
    Deletion,
    Insertion,
    Substitution,
    Variant,
    Variation,
    altbases,
    mutation,
    refbases,
    variations
using Statistics: mean
using VariantCallFormat: VCF
using XAM: BAM, SAM

export VariationCall
export VariationInfo
export VariationPileup
export altdepth
export call_variant
export depth
export filters
export frequency
export interval
export p_value
export pileup
export quality
export readpos
export relativepos
export seqpos
export strand
export strand_bias
export variation
export variation_test
export variationinfos
export vcf

include("fasta.jl")
include("xam.jl")
include("variation.jl")
include("interval.jl")
include("variationinfo.jl")
include("variationpileup.jl")
include("variationcall.jl")

const VERSION = project_version(
    string(joinpath(parent(parent(Path(Base.find_package("HapLink")))), "Project.toml"))
)

end #module
