module HapLink

using ArgParse
using BioAlignments
using BioSequences
using BioSymbols
using Combinatorics
using DataFrames
using Dates
using Distributions
using FASTX
using FilePaths
using HypothesisTests
using OrderedCollections
using Random
using SHA
using XAM
using YAML

export chromosome
export consensus_sequence
export consensus_variants
export frequency
export read_vcf

const VERSION = ArgParse.project_version(
    string(joinpath(parent(parent(Path(Base.find_package("HapLink")))), "Project.toml"))
)

include("variant.jl")
include("haplotype.jl")
include("readcounts.jl")
include("sequences.jl")
include("variantcalling.jl")
include("haplotypecalling.jl")
include("haplotypemeta.jl")
include("consensus.jl")

function parse_arguments()
    s = ArgParseSettings(;
        prog="haplink",
        description="A haplotype caller for long sequencing reads using linkage disequilibrium",
        version=VERSION,
        add_version=true,
        autofix_names=true,
    )
    # Disable Julia formatter as it doesn't understand the nested table syntax of ArgParse
    #! format: off

    # Declare the HapLink commands
    @add_arg_table s begin
        "variants"
            help = "Call variants"
            action = :command
        "consensus"
            help = "Create consensus sequence"
            action = :command
        "haplotypes"
            help = "Call haplotypes"
            action = :command
        "sequences"
            help = "Convert haplotypes to fasta"
            action = :command
    end #add_arg_table

    # Add arguments for the variant calling command
    @add_arg_table s["variants"] begin
        "--bam", "-i"
            help = """
                BAM formatted file containing the alignment to call haplotypes from
                """
            required     = true
            arg_type     = String
            range_tester = x -> isfile(x)
        "--reference", "-r"
            help = """
                FASTA formatted reference genome
                """
            required     = true
            arg_type     = String
            range_tester = x -> isfile(x)
        "--output", "-o"
            help = """
                File to output all variant calls to in VCF format
                """
            required     = true
        "--quality", "-q"
            help = """
                Minimum average quality (PHRED score) of a variant basecall
                """
            required     = false
            default      = 21
            arg_type     = Int64
            range_tester = x -> x > 0
        "--frequency", "-t"
            help = """
                Minimum frequency at which a variant must appear
                """
            required     = false
            default      = 0.05
            arg_type     = Float64
            range_tester = x -> (x >= 0) && (x <= 1)
        "--position", "-x"
            help = """
                Remove variants that occur only in positions within this percentage of the
                end
                """
            required     = false
            default      = 0.1
            arg_type     = Float64
            range_tester = x -> (x >= 0) && (x <= 1)
        "--significance", "-a"
            help = """
                Maximum Χ-squared significance level to consider a variant
                """
            required     = false
            default      = 1e-5
            arg_type     = Float64
            range_tester = x -> (x >= 0) && (x <= 1)
        "--depth", "-d"
            help = """
                Minimum depth to consider a variant
                """
            required     = false
            default      = 10
            arg_type     = Int64
            range_tester = x -> x >= 1
    end #add_arg_table

    # Add arguments for the consensus sequence generation command
    @add_arg_table! s["consensus"] begin
        "--reference", "-r"
            help = """
                FASTA formatted reference genome
                """
            required = true
            arg_type = String
            range_tester = x -> isfile(x)
        "--variants", "-v"
            help = """
                VCF formatted variant calls
                """
            required = true
            arg_type = String
            range_tester = x -> isfile(x)
        "--output", "-o"
            help = """
                File to output consensus sequence to in FASTA format
                """
            required = true
        "--frequency", "-t"
            help = """
                Minimum frequency at which a variant must appear to be included in the
                consensus sequence
                """
            required = false
            default = 0.5
            arg_type = Float64
            range_tester = x -> (x >= 0) && (x <= 1)
        "--prefix", "-p"
            help = "The prefix (sample identifier) to label the sequence as"
            required = false
            default = nothing
    end #add_arg_table

    # Add arguments for the haplotype calling command
    @add_arg_table s["haplotypes"] begin
        "--bam", "-i"
            help = """
                BAM formatted file containing the alignment to call haplotypes from
                """
            required     = true
            arg_type     = String
            range_tester = x -> isfile(x)
        "--variants", "-v"
            help = """
                VCF formmated file containing the variant calls to consider haplotypes from
                """
            required     = true
            arg_type     = String
            range_tester = x -> isfile(x)
        "--reference", "-r"
            help = """
                FASTA formatted reference genome
                """
            required     = true
            arg_type     = String
            range_tester = x -> isfile(x)
        "--output", "-o"
            help = """
                File to output all haplotype calls to in YAML format
                """
            required     = true
        "--significance", "-a"
            help = """
                Maximum Χ-squared significance level to consider a haplotype
                """
            required     = false
            default      = 1e-5
            arg_type     = Float64
            range_tester = x -> (x >= 0) && (x <= 1)
        "--depth", "-d"
            help = """
                Minimum depth to consider a haplotype
                """
            required     = false
            default      = 10
            arg_type     = Int64
            range_tester = x -> x >= 1
        "--method"
            help = """
                Haplotype read-building method. Choose 'ml-template' or 'raw'
                """
            required     = false
            default      = "ml-template"
            arg_type     = String
            range_tester = x -> (x == "ml-template") || (x == "raw")
        "--overlap_min", "-m"
            help = """
                The minimum amount reads must overlap to be placed in the same simulated
                template strand.
                """
            required = false
            default  = 0
            arg_type = Int64
        "--overlap_max", "-n"
            help = """
                The maximum amount reads are allowed to overlap to be placed in the same
                simulated template strand.
                """
            required     = false
            default      = 100
            arg_type     = Int64
            range_tester = x -> x >= 0
        "--iterations", "-j"
            help = """
                Formula to determine how many iterations to perform using one of the ml
                methods
                """
            required = false
            default  = 1000
            arg_type = Int64
        "--seed", "-s"
            help = """
                Set the random seed for the maximum likelihood psuedo-read generation
                process
                """
                required = false
                arg_type = UInt64
    end #add_arg_table

    @add_arg_table! s["sequences"] begin
        "--haplotypes", "-i"
            help = """
                YAML file describing the haplotypes to create
                """
                required     = true
                arg_type     = String
                range_tester = x -> isfile(x)
        "--reference", "-r"
            help = """
                FASTA formatted reference genome
                """
            required     = true
            arg_type     = String
            range_tester = x -> isfile(x)
        "--output", "-o"
            help = """
                File to output haplotype sequences to in FASTA format
                """
            required     = true
    end #add_arg_table
    #! format: on

    args = parse_args(s)

    return args["%COMMAND%"], args[args["%COMMAND%"]]
end #function

Base.@ccallable function haplink()::Cint
    command, arguments = parse_arguments()
    if command == "variants"
        variants(arguments)
    elseif command == "consensus"
        _consensus(arguments)
    elseif command == "haplotypes"
        haplotypes(arguments)
    elseif command == "sequences"
        sequences(arguments)
    else
        @error "Unknown command $command. Use 'variants', 'consensus', 'haplotypes', or 'sequences'."
        return 1
    end #if

    return 0
end #function

function variants(arguments::Dict{String,Any})
    # Read the argument table in as variables
    bamfile      = arguments["bam"]
    reffile      = arguments["reference"]
    outfile      = arguments["output"]
    quality      = arguments["quality"]
    frequency    = arguments["frequency"]
    position     = arguments["position"]
    significance = arguments["significance"]
    depth        = arguments["depth"]

    # Call variants
    variants = callvariants(
        countbasestats(bamfile, reffile), depth, quality, position, frequency, significance
    )

    # Save the variants to a VCF file
    return savevcf(variants, outfile, reffile, depth, quality, position, significance)
end #function

function _consensus(arguments::Dict{String,Any})
    reffile = Path(arguments["reference"])
    varfile = Path(arguments["variants"])
    outfile = Path(arguments["output"])
    frequency = arguments["frequency"]
    prefix = arguments["prefix"]

    if isnothing(prefix)
        prefix = filename(varfile)
    end #if

    variants = read_vcf(varfile)

    refseq = FASTA.sequence(_first_record(reffile))

    conseq = consensus_sequence(refseq, variants; freq=frequency)

    FASTA.Writer(open(outfile, "w")) do f
        write(f, FASTA.Record("$(prefix)_CONSENSUS", conseq))
    end #do

    return nothing
end #function

function haplotypes(arguments::Dict{String,Any})
    # Read the argument table in as variables
    bamfile = arguments["bam"]
    varfile = arguments["variants"]
    reffile = arguments["reference"]
    outfile = arguments["output"]
    significance = arguments["significance"]
    depth = arguments["depth"]
    method = arguments["method"]
    minoverlap = arguments["overlap_min"]
    maxoverlap = arguments["overlap_max"]
    iterations = arguments["iterations"]
    seed = arguments["seed"]

    # Set the random seed, if specified
    if !isnothing(seed)
        Random.seed!(seed)
    end #if

    # Read in the variants
    variants = read_vcf(varfile)

    if occursin("ml", method)
        # Create a read matching algorithm
        domatch =
            (r1::BAM.Record, r2::BAM.Record, pos::AbstractVecOrMat{Int}) ->
                BAM.position(r2) >= BAM.position(r1) &&
                    variant_positions_match(r1, r2, pos) &&
                    overlap_inrange(r1, r2; max=maxoverlap, min=minoverlap)

        # Use the simulated read method
        hapmethod =
            (h::Haplotype, b::AbstractString) ->
                simulate_genome(h, b; iterations=iterations, nextreadcandidates=domatch)
    else
        # Use the actual read method
        hapmethod = (h::Haplotype, b::AbstractString) -> longread_genome(h, b)
    end #if

    # Get the haplotypes and most of the statstics
    haplotypes = find_haplotypes(variants, bamfile, reffile, depth, significance, hapmethod)

    # Write the found haplotypes to file
    outdata = OrderedDict(
        "version" => VERSION,
        "settings" => arguments,
        "haplotypes" => Dict.(collect(haplotypes)),
    )
    YAML.write_file(outfile, outdata)

    return nothing
end #function

function sequences(arguments::Dict{String,Any})
    hfile = arguments["haplotypes"]
    rfile = arguments["reference"]
    ffile = arguments["output"]

    haplodata = YAML.load_file(hfile)["haplotypes"]
    consensus_snps = first(filter(d -> d["name"] == "CONSENSUS", haplodata))["snps"]
    consensus_hap = Haplotype(Variant.(consensus_snps))

    otherhaps = filter(d -> d["name"] != "CONSENSUS", haplodata)
    haplotypes = Haplotype.(map(f -> Variant.(f["snps"]), otherhaps))

    refrec = _first_record(rfile)

    conrec = FASTA.Record(
        FASTA.identifier(refrec),
        FASTA.description(refrec),
        mutate(FASTA.sequence(refrec), consensus_hap),
    )

    newrecords = unique(mutate.([conrec], haplotypes))

    open(FASTA.Writer, ffile) do f
        for r in newrecords
            write(f, r)
        end #for
    end #do

    return nothing
end #function

end #module
