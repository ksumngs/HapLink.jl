module HapLink

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args, project_version
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
    reconstruct!,
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

function _parse_arguments()
    s = ArgParseSettings(;
        prog="haplink",
        description="A haplotype caller for long sequencing reads using linkage disequilibrium",
        version=VERSION,
        add_version=true,
    )
    # Disable Julia formatter as it doesn't understand the nested table syntax of ArgParse
    #! format: off

    # Declare the HapLink commands
    @add_arg_table! s begin
        "variants"
            help = "Call variants"
            action = :command
        "consensus"
            help = "Create consensus sequence"
            action = :command
    end #add_arg_table

    # Add arguments for the variant calling command
    @add_arg_table! s["variants"] begin
        "reference"
            arg_type = String
            required = true
            range_tester = x -> isfile(x)
            help = "FASTA formatted reference genome sequence file"
        "bam"
            arg_type = String
            required = true
            range_tester = x -> isfile(x)
            help = "BAM formatted alignment file"
        "--output", "-o"
            arg_type = String
            help = "Write output to file (VCF format)"
        "--depth", "-m"
            arg_type = Int64
            default = 10
            range_tester = x -> x >= 1
            help = "Mimimum read depth to call variant"
        "--quality", "-q"
            arg_type = Float64
            default = 12.0
            range_tester = x -> x > 0
            help = "Minimum average basecall quality score to call variant"
        "--frequency", "-t"
            arg_type = Float64
            default = 0.05
            range_tester = x -> (x >= 0) && (x <= 1)
            help = "Minimum alternate base frequency to call variant"
        "--position", "-x"
            arg_type = Float64
            default = 0.1
            range_tester = x -> (x >= 0) && (x < 0.5)
            help = "Maximum distance from edge to call variant"
        "--strand", "-s"
            arg_type = Float64
            range_tester = x -> (x >= 0) && (x < 0.5)
            help = "Maximum proportion of alternate found on one strand to call variant"
        "--significance", "-a"
            arg_type = Float64
            default = 1e-5
            range_tester = x -> (x >= 0) && (x <= 1)
            help = "Maximum Χ-squared p-value level to call variant"
    end #add_arg_table

    @add_arg_table! s["consensus"] begin
        "reference"
            arg_type = String
            required = true
            range_tester = x -> isfile(x)
            help = "FASTA formatted reference genome sequence file"
        "variants"
            arg_type = String
            required = true
            range_tester = x -> isfile(x)
            help = "VCF formatted variant calls"
        "--output", "-o"
            arg_type = String
            help = "Write output to file (VCF format)"
        "--frequency", "-t"
            arg_type = Float64
            default = 0.5
            range_tester = x -> (x >= 0) && (x <= 1)
            help = "Minimum alternate base frequency to label variant as consensus"
        "--prefix", "-p"
            arg_type = String
            help = "The prefix (sequence identifier) of the consensus sequence"
    end #add_arg_table
    #! format: on

    args = parse_args(s)

    return args["%COMMAND%"], args[args["%COMMAND%"]]
end #function

Base.@ccallable function julia_main()::Cint
    cmd, args = _parse_arguments()

    if cmd == "variants"
        _haplink_variants(args)
    elseif cmd == "consensus"
        _haplink_consensus(args)
    else
        return 1
    end #if

    return 0
end #function

function _haplink_variants(args::Dict{String,Any})
    reffile = args["reference"]
    bam = args["bam"]
    outfile = args["output"]
    depth = args["depth"]
    qual = args["quality"]
    freq = args["frequency"]
    pos = args["position"]
    strand = args["strand"]
    α = args["significance"]

    refname = FASTA.identifier(_first_record(reffile))

    vcf_head = _vcf_header(reffile, α; D=depth, Q=qual, X=pos, F=freq, S=strand)
    out_stream = isnothing(outfile) ? stdout : open(outfile, "w")
    vcf_writer = VCF.Writer(out_stream, vcf_head)

    pile = pileup(bam, reffile)

    for p in pile
        vc = call_variant(p, α; D=depth, Q=qual, X=pos, F=freq, S=strand)
        vcf_rec = vcf(vc, refname)
        write(vcf_writer, vcf_rec)
    end #for

    close(out_stream)

    return 0
end #function

function _haplink_consensus(args::Dict{String,Any})
    reffile = args["reference"]
    varfile = args["variants"]
    outfile = args["output"]
    freq = args["frequency"]
    prefix = args["prefix"]

    refrec = _first_record(reffile)
    ref_id = isnothing(prefix) ? FASTA.identifier(refrec) : prefix
    ref_seq = FASTA.sequence(refrec)

    vars = Variation{typeof(ref_seq),eltype(ref_seq)}[]

    vcf_reader = VCF.Reader(open(varfile, "r"))
    for vcf_rec in vcf_reader
        parse(Float64, VCF.info(vcf_rec, "AF")) >= freq || continue
        first(VCF.filter(vcf_rec)) == "PASS" || continue
        push!(vars, variation(vcf_rec, ref_seq))
    end #for
    close(vcf_reader)

    con_seq = isempty(vars) ? ref_seq : reconstruct!(ref_seq, Variant(ref_seq, vars))

    FASTA.Writer(isnothing(outfile) ? stdout : open(outfile, "w")) do fasta_writer
        fasta_record = FASTA.Record("$(ref_id)_CONSENSUS", con_seq)
        write(fasta_writer, fasta_record)
    end #do

    return 0
end #function

end #module
