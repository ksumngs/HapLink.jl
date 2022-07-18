module HapLink

using ArgParse: ArgParseSettings, @add_arg_table!, project_version
using BioAlignments: Alignment, AlignedSequence, PairwiseAlignment, ref2seq
using BioGenerics: BioGenerics, leftposition, rightposition, metadata
using BioSequences: BioSequence, NucleotideSeq
using BioSymbols: BioSymbol
using FASTX: FASTA
using FilePaths: FilePaths, AbstractPath, Path
using GenomicFeatures: Interval, Strand, STRAND_POS, eachoverlap
using SequenceVariation:
    SequenceVariation,
    Deletion,
    Insertion,
    Substitution,
    Variant,
    Variation,
    mutation,
    variations
using Statistics: mean
using XAM: BAM, SAM

export VariationCall
export VariationInfo
export VariationPileup
export altdepth
export depth
export frequency
export interval
export pileup
export quality
export readpos
export relativepos
export seqpos
export strand
export strand_bias
export variation
export variationinfos

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
