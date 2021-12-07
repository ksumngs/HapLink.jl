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

export countbasestats
export callvariants
export findsimulatedoccurrences
export linkage
export sumsliced
export mutate

include("variant.jl")
include("haplotype.jl")
include("readcounts.jl")
include("sequences.jl")
include("variantcalling.jl")
include("haplotypecalling.jl")

Base.@ccallable function haplink()::Cint
    s = ArgParseSettings(
        prog="haplink",
        description="A haplotype caller for long sequencing reads using linkage disequilibrium",
        version=VERSION,
        add_version=true,
        autofix_names=true
    )
    @add_arg_table s begin
        "bamfile"
            help         = "BAM formatted file containing the alignment to call haplotypes from"
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
            help         = "(NOT IMPLEMENTED) File to output all variant calls in VCF format"
            required     = false
        "--prefix", "-p"
            help         = "Test to start the output file names with. If unspecified, will use the name of the alignment file up to the first dot (.)"
            required     = false
            arg_type     = String
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
            help         = "Haplotype read-building method. Choose one of 'ml-overlap', 'ml-gapped' or 'raw'"
            required     = false
            default      = "ml-overlap"
            arg_type     = String
            range_tester = x -> (x == "ml-overlap") || (x == "ml-gapped") || (x == "raw")
        "--iterations"
            help         = "Formula to determine how many iterations to perform using one of the ml methods"
            required     = false
            default      = ":(1000)"
            arg_type     = String
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
        α_variant
    )

    # Save the variants to a VCF file, if requested
    if !isnothing(args["variants"])
        savevcf(
            variants,
            args["variants"],
            reffile,
            D_variant,
            Q_variant,
            x_variant,
            α_variant
        )
    end #if

    # Check for zero found variants
    if length(variants) < 1
        @warn "No variants found!"
        touch(string(prefix, ".yaml"))
        cp(reffile, string(prefix, ".fasta"))
        return
    end #if

    if occursin("ml", args["method"])
        # TODO: implement an expression-evaluator for ML iterations
        # Calculate the number of iterations for each haplotype
        iterations = 1000 # max(1000, D_haplotype*length(variants)^2)

        # TODO: implement an overlapped-read ML algorithm

        haplotypes = findsimulatedhaplotypes(
            variants,
            bamfile,
            D_haplotype,
            α_haplotype,
            iterations=iterations
        )
    else
        # TODO: implement a long read, non-simulated likage finder
        haplotypes = findhaplotypes(variants, bamfile, D_haplotype, α_haplotype)
    end #if

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

"""
    function findsimulatedoccurrences(args...; kwargs...)

Find the number of times a particular haplotype is supported by a maximum likelihood
simulation of combining aligned reads in a BAM file.

# Arguments
- `haplotype::Haplotype`: The combination of variants to test the aligned reads for evidence
    of
- `bamfile::AbstractString`: The path to a BAM file containing aligned reads to be tested
    for instances of `haplotype`

# Keywords
- `iterations::Integer=1000`: The number of times to combine reads and test for the presence
    of `haplotype`

`findsimulatedoccurrences` examines each variant position in `haplotype` and finds a random
read containing that position from `bamfile`. It will then check if the next variant
position is contained on the previous read, and if not pick a new random read that contains
that variant position. In this way, it assembles a set of reads that conceivably could have
come from the same template strand via maximum likelihood.

From that set of reads, `findsimulatedoccurrences` returns an ``N``-dimensional matrix where
``N`` is the number of variant positions in `haplotypes`. The ``1`` index position in the
``n``th dimension represents the number of times the ``n``th variant position was found to
have the reference base called, while the ``2`` index position represents the number of
times the ``n``th variant position was found to have the alternate base called. E.g.
`first(findsimulatedoccurrences(...))` gives the number of times the all-reference base
haplotype was found in the simulation, while `findsimulatedoccurrences(...)[end]` gives the
number of times the all-alternate base haplotype was found.
"""
function findsimulatedoccurrences(
    haplotype::Haplotype,
    bamfile::AbstractString;
    iterations=1000
)

    # Extract the SNPs we care about
    mutations = haplotype.mutations

    # Create an empty array for the simulated long reads
    pseudoreads = Array{Symbol}(undef, iterations, length(mutations))

    # Start reading the BAM file
    open(BAM.Reader, bamfile) do bamreader
        # Collect the reads
        reads = collect(bamreader)

        # Start iterating
        Threads.@threads for i ∈ 1:iterations
            # Get the reads that contain the first mutation
            lastcontainingreads = filter(
                b -> BAM.position(b) < mutations[1].position && BAM.rightposition(b) > mutations[1].position,
                reads
            )

            # Pull a random read from that pool
            lastread = rand(lastcontainingreads)

            # Find this read's basecall at that position
            basecall = baseatreferenceposition(lastread, mutations[1].position)
            basematch = matchvariant(basecall, mutations[1])

            pseudoreads[i, 1] = basematch

            for j ∈ 2:length(mutations)
                if (BAM.position(lastread) < mutations[j].position && BAM.rightposition(lastread) > mutations[j].position)
                    thisread = lastread
                else
                    thiscontainingreads = filter(
                        b -> BAM.position(b) > BAM.rightposition(lastread) && BAM.position(b) < mutations[j].position && BAM.rightposition(b) > mutations[j].position,
                        reads
                    )
                    if length(thiscontainingreads) < 1
                        pseudoreads[i,j] = :other
                        continue
                    end #if
                    thisread = rand(thiscontainingreads)
                end #if

                # Find this read's basecall at that position
                basecall = baseatreferenceposition(thisread, mutations[j].position)
                basematch = matchvariant(basecall, mutations[j])

                pseudoreads[i, j] = basematch

                lastread = thisread
            end #for
        end #for
    end #do

    # Set up haplotype counts
    hapcounts = zeros(Int, repeat([2], length(mutations))...)

    for i ∈ 1:iterations
        matches = pseudoreads[i, :]
        if !any(matches .== :other)
            coordinate = CartesianIndex((Int.(matches .== :alternate) .+ 1)...)
            hapcounts[coordinate] += 1
        end #if
    end #for

    return hapcounts
end #function

"""
    linkage(counts::AbstractArray{Int})

Calculates the linkage disequilibrium and Chi-squared significance level of a combination of
haplotypes whose number of occurrences are given by `counts`.

`counts` is an ``N``-dimensional array where the ``N``th dimension represents the ``N``th
variant call position within a haplotype. `findoccurrences` produces such an array.
"""
function linkage(counts::AbstractArray{Int})
    # Get the probability of finding a perfect reference sequence
    P_allref = first(counts) / sum(counts)

    # Get the probabilities of finding reference bases in any of the haplotypes
    P_refs = sumsliced.([counts], 1:ndims(counts)) ./ sum(counts)

    # Calculate linkage disequilibrium
    Δ = P_allref - prod(P_refs)

    # Calculate the test statistic
    r = Δ / (prod(P_refs .* (1 .- P_refs))^(1/ndims(counts)))
    Χ_squared = r^2 * sum(counts)

    # Calculate the significance
    p = 1 - cdf(Chisq(1), Χ_squared)

    return Δ, p
end #function

"""
    sumsliced(A::AbstractArray, dim::Int, pos::Int=1)

Sum all elements that are that can be referenced by `pos` in the `dim` dimension of `A`.

# Example

```jldoctest
julia> A = reshape(1:8, 2, 2, 2)
2×2×2 reshape(::UnitRange{Int64}, 2, 2, 2) with eltype Int64:
[:, :, 1] =
 1  3
 2  4

[:, :, 2] =
 5  7
 6  8

julia> sumsliced(A, 2)
14

julia> sumsliced(A, 2, 2)
22
```

Heavily inspired by Holy, Tim "Multidimensional algorithms and iteration"
<https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes>
"""
function sumsliced(A::AbstractArray, dim::Int, pos::Int=1)
    i_pre  = CartesianIndices(size(A)[1:dim-1])
    i_post = CartesianIndices(size(A)[dim+1:end])
    return sum(A[i_pre, pos, i_post])
end #function

function serialize_yaml(h::Pair; reason::Union{String,Nothing}=nothing)
    occurrences = "occurrences:\n"
    for i in CartesianIndices(h.second)
        location = [Tuple(i)...]
        variantpattern = string.(replace(replace(location, 1 => "ref"), 2 => "alt"))
        key = join(variantpattern, "_")
        occurrences = string(occurrences, "  ", key, ": ", h.second[i], "\n")
    end #for

    return string(
        serialize_yaml(h.first, reason=reason),
        occurrences,
        "Δ: ",
        linkage(h.second)[1],
        "\n",
        "p: ",
        linkage(h.second)[2],
        "\n"
    )
end #function

end #module
