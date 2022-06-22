module HapLink

using ArgParse: ArgParseSettings, @add_arg_table!, project_version
using SequenceVariation:
    SequenceVariation,
    Deletion,
    Insertion,
    Substitution,
    Variant,
    Variation,
    mutation,
    variations
using BioSymbols: BioSymbol
using XAM: BAM, SAM
using BioAlignments: Alignment, AlignedSequence, PairwiseAlignment, ref2seq
using BioGenerics: BioGenerics, leftposition
using FilePaths: Path

export relativepos
export seqpos

_xam_switch(r::Union{Type{SAM.Record},Type{BAM.Record}}) = r <: SAM.Record ? :SAM : :BAM

include("variation.jl")

const VERSION = project_version(
    string(joinpath(parent(parent(Path(Base.find_package("HapLink")))), "Project.toml"))
)

end #module
