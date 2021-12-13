module HapLink

using ArgParse
using FASTX
using Dates
using HypothesisTests
using DataFrames
using FilePaths
using BioSequences
using BioAlignments
using BioSymbols
using Combinatorics
using Distributions
using SHA
using XAM
using YAML

const VERSION = ArgParse.project_version(
    string(joinpath(parent(parent(Path(Base.find_package("HapLink")))), "Project.toml"))
)

include("variant.jl")
include("haplotype.jl")
include("readcounts.jl")
include("sequences.jl")
include("variantcalling.jl")
include("haplotypecalling.jl")

Base.@ccallable function haplink()::Cint
    s = ArgParseSettings(;
        prog="haplink",
        description="A haplotype caller for long sequencing reads using linkage disequilibrium",
        version=VERSION,
        add_version=true,
        autofix_names=true,
    )
    @add_arg_table s begin
        "bamfile"
        help = "BAM formatted file containing the alignment to call haplotypes from"
        required     = true
        arg_type     = String
        range_tester = x -> isfile(x)
        "--reference", "-r"
        help         = "FASTA formatted reference genome"
        required     = true
        arg_type     = String
        range_tester = x -> isfile(x)
        "--annotations", "-g"
        help         = "GFF3 formatted annotations for the reference genome"
        required     = false
        arg_type     = String
        range_tester = x -> isfile(x)
        "--variants", "-v"
        help     = "(NOT IMPLEMENTED) File to output all variant calls in VCF format"
        required = false
        "--prefix", "-p"
        help     = "Test to start the output file names with. If unspecified, will use the name of the alignment file up to the first dot (.)"
        required = false
        arg_type = String
        "--quality", "-q"
        help         = "Minimum average quality (PHRED score) of a variant basecall"
        required     = false
        default      = 21
        arg_type     = Int64
        range_tester = x -> x > 0
        "--frequency", "-t"
        help         = "Minimum frequency at which a variant must appear"
        required     = false
        default      = 0.05
        arg_type     = Float64
        range_tester = x -> (x >= 0) && (x <= 1)
        "--position", "-x"
        help         = "Remove variants that occur only in positions within this percentage of the end"
        required     = false
        default      = 0.1
        arg_type     = Float64
        range_tester = x -> (x >= 0) && (x <= 1)
        "--variant_significance", "-a"
        help         = "Maximum Χ-squared significance level to consider a haplotype"
        required     = false
        default      = 1e-5
        arg_type     = Float64
        range_tester = x -> (x >= 0) && (x <= 1)
        "--haplotype_significance", "-b"
        help         = "Maximum Χ-squared significance level to consider a haplotype"
        required     = false
        default      = 1e-5
        arg_type     = Float64
        range_tester = x -> (x >= 0) && (x <= 1)
        "--variant_depth", "-d"
        help         = "Minimum depth to consider a variant"
        required     = false
        default      = 10
        arg_type     = Int64
        range_tester = x -> x >= 1
        "--haplotype_depth", "-u"
        help         = "Minimum depth to consider a haplotype"
        required     = false
        default      = 10
        arg_type     = Int64
        range_tester = x -> x >= 1
        "--method"
        help = "Haplotype read-building method. Choose one of 'ml-overlap', 'ml-gapped' or 'raw'"
        required     = false
        default      = "ml-overlap"
        arg_type     = String
        range_tester = x -> (x == "ml-overlap") || (x == "ml-gapped") || (x == "raw")
        "--overlap_min", "-i"
        help     = """
                The minimum amount reads must overlap to be placed in the same simulated \
                template strand.
                """
        required = false
        default  = 0
        arg_type = Int64
        "--overlap_max", "-m"
        help     = """
                The maximum amount reads are allowed to overlap to be placed in the same \
                simulated template strand.
                """
        required = false
        default  = 100
        arg_type = Int64
        range_tester = x -> x >= 0
        "--iterations"
        help = "Formula to determine how many iterations to perform using one of the ml methods"
        required = false
        default  = ":(1000)"
        arg_type = String
    end

    args = parse_args(s)

    # 1. Analyze bam
    # 2. Call variants
    # 3. Call haplotypes
    # 4. Export haplotypes as YAML
    # 5. Export haplotypes as FASTA

    # Read the argument table in as variables
    bamfile        = args["bamfile"]
    reffile        = args["reference"]
    annotationfile = args["annotations"]
    Q_variant      = args["quality"]
    f_variant      = args["frequency"]
    x_variant      = args["position"]
    α_variant      = args["variant_significance"]
    α_haplotype    = args["haplotype_significance"]
    D_variant      = args["variant_depth"]
    D_haplotype    = args["haplotype_depth"]
    maxoverlap     = args["overlap_max"]
    minoverlap     = args["overlap_min"]

    # Find the file prefix for output files if none was provided
    bampath = Path(bamfile)
    prefix = isnothing(args["prefix"]) ? filename(bampath) : args["prefix"]

    # Call variants
    variants = callvariants(
        countbasestats(bamfile, reffile),
        D_variant,
        Q_variant,
        x_variant,
        f_variant,
        α_variant,
    )

    # Save the variants to a VCF file, if requested
    if !isnothing(args["variants"])
        savevcf(
            variants, args["variants"], reffile, D_variant, Q_variant, x_variant, α_variant
        )
    end #if

    # Check for zero found variants
    if length(variants) < 1
        @warn "No variants found!"
        touch(string(prefix, ".yaml"))
        cp(reffile, string(prefix, ".fasta"))
        return 0
    end #if

    if occursin("ml", args["method"])
        # Calculate the number of iterations for each haplotype
        # TODO: implement an expression-evaluator for ML iterations
        iterations = 1000

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

    haplotypes = find_haplotypes(variants, bamfile, D_haplotype, α_haplotype, hapmethod)

    # Write the found haplotypes to file
    yamlfile = string(prefix, ".yaml")
    open(string(prefix, ".yaml"), "w") do f
        for happair in haplotypes
            write(f, serialize_yaml(happair))
        end #for
    end #do

    open(FASTA.Reader, reffile) do r
        record = collect(r)[1]
        newrecords = unique(mutate.([record], collect(keys(haplotypes))))

        # Write the found haplotypes to FASTA
        open(FASTA.Writer, string(prefix, ".fasta")) do f
            for q in newrecords
                write(f, q)
            end #for
        end #do
    end #do

    return 0
end #function

Base.@ccallable function make_haplotype_fastas()::Cint
    hfile = ARGS[1]
    rfile = ARGS[2]
    ffile = ARGS[3]

    haployaml = read(hfile, String)
    haplostrings = split(haployaml, "---\n")[2:end]
    haplotypes = Haplotype.(YAML.load.(haplostrings, dicttype=Dict{String,Any}))

    rreader = open(FASTA.Reader, rfile)
    refrec = collect(rreader)[1]
    close(rreader)

    newrecords = unique(mutate.([refrec], haplotypes))

    fwriter = open(FASTA.Writer, ffile)
    for r in newrecords
        write(fwriter, r)
    end
    close(fwriter)

    return 0
end #function

end #module
