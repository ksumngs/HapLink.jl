module HapLink

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args, project_version
using BioAlignments: Alignment, AlignedSequence, PairwiseAlignment, ref2seq
using BioGenerics: BioGenerics, leftposition, rightposition, metadata
using BioSequences: BioSequence, LongDNA, NucleotideSeq
using BioSymbols: BioSymbol
using Combinatorics: combinations
using Dates: Dates, today
using Distributions: Chisq, cdf
using FASTX: FASTA
using FilePaths: FilePaths, AbstractPath, Path, absolute
using GenomicFeatures: Interval, Strand, STRAND_POS, eachoverlap
using HypothesisTests: FisherExactTest, pvalue
using OrderedCollections: OrderedDict
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
    reference,
    variations
using SHA: sha1
using Statistics: mean
using VariantCallFormat: VCF
using XAM: BAM, SAM
using YAML: YAML

export VariationCall
export VariationInfo
export VariationPileup
export altdepth
export call_variant
export consensus
export depth
export filters
export findset
export frequency
export interval
export linkage
export occurence_matrix
export p_value
export pileup
export quality
export readpos
export relativepos
export seqpos
export strand
export strand_bias
export subconsensus_variations
export sumsliced
export variation
export variation_test
export variationinfos
export vcf

include("fasta.jl")
include("xam.jl")
include("variant.jl")
include("variation.jl")
include("interval.jl")
include("variationinfo.jl")
include("variationpileup.jl")
include("variationcall.jl")
include("consensus.jl")
include("findset.jl")
include("haplotypecalling.jl")

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
        "haplotypes"
            help = "Call haplotypes"
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

    @add_arg_table! s["haplotypes"] begin
        "reference"
            arg_type = String
            required = true
            range_tester = isfile
            help = "FASTA formatted reference genome sequence file"
        "variants"
            arg_type = String
            required = true
            range_tester = isfile
            help = "VCF formatted variant calls"
        "bam"
            arg_type = String
            required = true
            range_tester = isfile
            help = "BAM formatted alignment file"
        "--output", "-o"
            arg_type = String
            help = "Write output to file (YAML format)"
        "--consensus_frequency", "-c"
            arg_type = Float64
            default = 0.5
            range_tester = x -> (x >= 0) && (x <= 1)
            help = "Minimum alternate base frequency to label variant as consensus"
        "--significance", "-a"
            arg_type = Float64
            default = 1e-5
            range_tester = x -> (x >= 0) && (x <= 1)
            help = "Maximum Χ-squared significance to call a haplotype"
        "--depth", "-d"
            arg_type = Int64
            default = 10
            range_tester = x -> x >= 1
            help = "Minimum depth to call a haplotype"
        "--frequency", "-t"
            arg_type = Float64
            default = 0.1
            range_tester = x -> (x >= 0) && (x <= 1)
            help = "Minimum frequency to call a haplotype"
        "--full_length_reads"
            action = :store_true
            help = "Use reads that span the entire genome to call haplotypes"
        "--simulated_reads"
            action = :store_true
            help = "Simulate long reads based on overlapping alignments. Overrides `--full_length_reads`"
        "--overlap_min", "-m"
            arg_type = Int64
            default = 0
            help = "The minimum number of bases that aligned reads must overlap to be compiled into the same simulated read. A value of `0` indicates that all gaps are permitted, while a negative value puts a cap on the length of gaps permitted. Only applies if `--simulated_reads` is passed."
        "--overlap_max", "-n"
            arg_type = Int64
            default = 500
            help = "The maximum number of bases that aligned reads are allowed to overlap to be compiled into the same simulated read. Only applies if `--simulated_reads` is passed."
        "--iterations", "-j"
            arg_type = UInt64
            default = UInt64(10000)
            help = "Number of simulated reads to attepmt to generate. Only applies if `--simulated_reads` is passed."
        "--seed", "-s"
            arg_type = UInt64
            help = "Set the random seed for the simulated read generation process. Only applies if `--simulated_reads` is passed."
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
    elseif cmd == "haplotypes"
        _haplink_haplotypes(args)
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

    consensus_record = consensus(reffile, varfile; frequency=freq, prefix=prefix)

    FASTA.Writer(isnothing(outfile) ? stdout : open(outfile, "w")) do fasta_writer
        write(fasta_writer, consensus_record)
    end #do

    return 0
end #function

function _haplink_haplotypes(args::Dict{String,Any})
    reffile = args["reference"]
    varfile = args["variants"]
    bamfile = args["bam"]
    outfile = args["output"]
    frequency = args["frequency"]
    consensus_frequency = args["consensus_frequency"]
    significance = args["significance"]
    depth = args["depth"]
    simulate_reads = args["simulated_reads"]

    refrecord = _first_record(reffile)
    refseq = FASTA.sequence(LongDNA{4}, refrecord)

    consensus_variant = consensus(refseq, varfile; frequency=consensus_frequency)

    consensus_sequence = reconstruct!(refseq, consensus_variant)

    read_pool = Pseudoread[]
    if simulate_reads
        @warn "Not implemented yet"
    else
        read_pool = _read.(_pseudoreads(bamfile, consensus_sequence))
    end #if

    subconsensus_vars = subconsensus_variations(varfile, consensus_variant)

    haplotype_tester =
        h -> ishaplotype(
            h,
            read_pool;
            min_frequency=frequency,
            significance_level=significance,
            min_depth=depth,
        )

    valid_haplotypes = findset(subconsensus_vars, haplotype_tester)

    hapcalls = HaplotypeCall[]

    for hap in valid_haplotypes
        push!(hapcalls, HaplotypeCall(hap, read_pool))
    end #for

    outdict = OrderedDict{String,Any}()
    outdict["version"] = VERSION
    outdict["settings"] = args
    outdict["haplotypes"] = [_dict(h) for h in hapcalls]

    out_stream = isnothing(outfile) ? stdout : open(outfile, "w")
    YAML.write(out_stream, outdict)

    return 0
end #function

end #module
