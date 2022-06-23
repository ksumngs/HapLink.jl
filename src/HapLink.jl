module HapLink

using ArgParse: ArgParseSettings, @add_arg_table!, project_version
using BioAlignments: Alignment, AlignedSequence, PairwiseAlignment, ref2seq
using BioGenerics: BioGenerics, leftposition, rightposition
using BioSequences: BioSequence, NucleotideSeq
using BioSymbols: BioSymbol
using FilePaths: Path
using GenomicFeatures: Interval, Strand
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

export VariationInfo
export quality
export readpos
export relativepos
export seqpos
export strand
export variation
export variationinfos

include("xam.jl")
include("variation.jl")
include("variationinfo.jl")

const VERSION = project_version(
    string(joinpath(parent(parent(Path(Base.find_package("HapLink")))), "Project.toml"))
)

end #module
