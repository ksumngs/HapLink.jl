module HapLink

using ArgParse: ArgParseSettings, @add_arg_table!, project_version
using BioAlignments: Alignment, AlignedSequence, PairwiseAlignment, ref2seq
using BioGenerics: BioGenerics, leftposition, rightposition
using BioSequences: BioSequence, NucleotideSeq
using BioSymbols: BioSymbol
using FilePaths: Path
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

export quality
export relativepos
export seqpos

_xam_switch(r::Union{Type{SAM.Record},Type{BAM.Record}}) = r <: SAM.Record ? :SAM : :BAM

include("variation.jl")

const VERSION = project_version(
    string(joinpath(parent(parent(Path(Base.find_package("HapLink")))), "Project.toml"))
)

end #module
