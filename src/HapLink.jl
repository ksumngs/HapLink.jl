"""
    HapLink

Viral haplotype calling via linkage disequilibrium
"""
module HapLink

using BioAlignments: Alignment, AlignedSequence, PairwiseAlignment, ref2seq
using BioGenerics: BioGenerics, leftposition, rightposition, metadata
using BioSequences: BioSequence, LongDNA, NucleotideSeq
using BioSymbols: BioSymbol, DNA
using Combinatorics: combinations
using Comonicon: Arg, @cast, @main, get_version
using Dates: Dates, today
using Distributions: Chisq, cdf
using FASTX: FASTA
using FilePaths: FilePaths, AbstractPath, Path, absolute
using GenomicFeatures: Interval, Strand, STRAND_POS, eachoverlap
using HypothesisTests: FisherExactTest, pvalue
using OrderedCollections: OrderedDict
using Random: seed!
using SequenceVariation:
    SequenceVariation,
    Deletion,
    Haplotype,
    Insertion,
    Substitution,
    Variation,
    altbases,
    mutation,
    reconstruct,
    refbases,
    reference,
    translate,
    variations
using SHA: sha1
using Statistics: mean
using VariantCallFormat: VCF
using XAM: BAM, SAM
using YAML: YAML

export Pseudoread
export VariationCall
export VariationInfo
export VariationPileup
export altdepth
export call_variant
export cigar
export consensus_haplotype
export consensus_record
export contridicts
export depth
export filters
export findset
export frequency
export interval
export linkage
export magnitude
export occurence_matrix
export overlap
export overlap_inrange
export overlapping_variations
export p_value
export pileup
export pseudoreads
export quality
export readpos
export relativepos
export seqpos
export significance
export simulate
export strand
export strand_bias
export subconsensus_variations
export sumsliced
export variant
export variation
export variation_test
export variationinfos
export variations_match
export vcf

# Declare some handy math aliases
const Σ = sum
const Π = prod

# Get the project version
const HAPLINK_VERSION = get_version(HapLink)

include("fasta.jl")
include("xam.jl")
include("haplotype.jl")
include("variation.jl")
include("interval.jl")
include("variationinfo.jl")
include("variationpileup.jl")
include("variationcall.jl")
include("consensus.jl")
include("pseudoread.jl")
include("findset.jl")
include("haplotypecalling.jl")

"""
    function variants(
        reference::String,
        bam::String;
        outfile::Union{String,Nothing}=nothing,
        significance::Float64=1e-5,
        depth::UInt64=UInt64(10),
        quality::Float64=12.0,
        frequency::Float64=0.05,
        position::Float64=0.5,
        strandedness::Union{Float64,Nothing}=nothing,
    )

Call variants

# Introduction

Decides which variations found within an alignment are real, and which are due to random
chance. HapLink uses Fisher's Exact Test to determine the statistical significance of
sequence variations, and optionally allows for other thresholds to reduce random noise in
the variant calling. Outputs a Variant Call Format (VCF) file compliant with VCF v4.

# Arguments

- `reference`: path to the reference genome to call variants against in fasta format. Must
    not be gzipped, but does not need to be indexed (have a sidecar fai file). HapLink only
    supports single-segment reference genomes: if `reference` includes more than one
    sequence, all but the first will be ignored.
- `bam`: alignment file to call variants from. Can be in SAM or BAM format, and does not
    need to be sorted or indexed, but variant calling speed will increase significantly if
    using a sorted and indexed (has a sidebar bai file) BAM file.

# Options

- `--outfile=<path>`: The file to write variant calls to. If left blank, variant calls are
    written to standard output.
- `--significance=<float>`: The alpha value for statistical significance of variant calls.
- `--depth=<int>`: Minimum number of times the variation must be observed within the
    alignment to be called a variant
- `--quality=<float>`: The minimum average basecall quality score for a variation to be
    called a variant
- `--frequency=<float>`: The minimum proportion of reads that the variation must be observed
    within compared to all reads covering its position for that variation to be called a
    variant
- `--position=<float>`: The distance (as a percentage) from the edge of reads that a
    variation must be observed at to be called a variant
- `--strandedness=<float>`: The maximum proportion of times that a variation can be observed
    on one strand versus the other to be called a variant. This metric is totally useless on
    single-stranded sequencing protocols like Oxford Nanopore, but can be useful for
    combining data between stranded protocols like most Illumina and Pacific Bio.
"""
@cast function variants(
    reference::String,
    bam::String;
    outfile::Union{String,Nothing}=nothing,
    significance::Float64=1e-5,
    depth::UInt64=UInt64(10),
    quality::Float64=12.0,
    frequency::Float64=0.05,
    position::Float64=0.5,
    strandedness::Union{Float64,Nothing}=nothing,
)
    refname = FASTA.identifier(_first_record(reference))

    vcf_head = _vcf_header(
        reference, significance; D=depth, Q=quality, X=position, F=frequency, S=strandedness
    )
    out_stream = isnothing(outfile) ? stdout : open(outfile, "w")
    vcf_writer = VCF.Writer(out_stream, vcf_head)

    pile = pileup(bam, reference)

    for p in pile
        vc = call_variant(
            p, significance; D=depth, Q=quality, X=position, F=frequency, S=strandedness
        )
        vcf_rec = vcf(vc, refname)
        write(vcf_writer, vcf_rec)
    end #for

    close(out_stream)

    return 0
end #function

"""
    function consensus(
        reference::String,
        variants::String;
        outfile::Union{String,Nothing}=nothing,
        frequency::Float64=0.5,
        prefix::String=""
    )

Convert variant calls to consensus sequence

# Introduction

Generates a consensus sequence based on a reference genome and previously called variants.
Will only consider variants with a "PASS" filter to be able to contribute to the consensus.
Outputs results in FASTA format.

# Arguments

- `reference`: The path to the reference genome. Must be in FASTA format.
- `variants`: The path to the variant calls. Must be in VCF v4 format.

# Options

- `--outfile=<path>`: The file to write the consensus sequence to. If left blank, the
    consensus sequence is written to standard output.
- `--frequency=<float>`: The minimum frequency at which a variant must appear to be
    considered part of the consensus. Note that HapLink does not support ambigous base
    calling (e.g. `N`, `R`, `Y`, etc.) at low frequencies unlike many variant callers.
- `--prefix=<string>`: Name of the new sequence. Defaults to using the FASTA identifier of
    the reference sequence.
"""
@cast function consensus(
    reference::String,
    variants::String;
    outfile::Union{String,Nothing}=nothing,
    frequency::Float64=0.5,
    prefix::String="",
)
    FASTA.Writer(isnothing(outfile) ? stdout : open(outfile, "w")) do fasta_writer
        write(
            fasta_writer,
            consensus_record(reference, variants; frequency=frequency, prefix=prefix),
        )
    end #do

    return 0
end #function

function _haplink_haplotypes(args::Dict{String,Any})
    reffile = args["reference"]
    varfile = args["variants"]
    bamfile = args["bam"]
    outfile = args["output"]
    consensus_frequency = args["consensus_frequency"]
    significance = args["significance"]
    depth = args["depth"]
    frequency = args["frequency"]
    simulate_reads = args["simulated_reads"]
    overlap_min = args["overlap_min"]
    overlap_max = args["overlap_max"]
    iterations = args["iterations"]
    input_seed = args["seed"]

    refrecord = _first_record(reffile)
    refseq = FASTA.sequence(LongDNA{4}, refrecord)

    consensus_variant = consensus_haplotype(refseq, varfile; frequency=consensus_frequency)
    consensus_sequence = reconstruct(consensus_variant)
    consensus_alignment = PairwiseAlignment(
        AlignedSequence(consensus_sequence, Alignment(cigar(consensus_variant))), refseq
    )

    fake_reads = pseudoreads(bamfile, consensus_sequence)

    subconsensus_vars = subconsensus_variations(varfile, consensus_variant)
    subconsensus_remapped = map(v -> translate(v, consensus_alignment), subconsensus_vars)
    sort!(subconsensus_remapped)

    read_pool = Haplotype{LongDNA{4},DNA}[]
    if simulate_reads
        if !isnothing(input_seed)
            seed!(input_seed)
        end #if

        small_pool = Vector{Union{Haplotype{LongDNA{4},DNA},Missing}}(undef, iterations)

        Threads.@threads for i in 1:iterations
            small_pool[i] = simulate(
                fake_reads,
                consensus_sequence,
                subconsensus_remapped;
                reverse_order=iseven(i),
                overlap_min=overlap_min,
                overlap_max=overlap_max,
            )
        end #for

        read_pool = collect(skipmissing(small_pool))
    else
        first_pos = leftposition(first(subconsensus_remapped))
        last_pos = rightposition(last(subconsensus_remapped))

        for fr in fake_reads
            leftposition(fr) >= first_pos || continue
            rightposition(fr) <= last_pos || continue

            v = variant(fr)
            t = translate(v, consensus_alignment)
            push!(read_pool, t)
        end #for
    end #if

    haplotype_tester =
        h -> ishaplotype(
            h,
            read_pool;
            min_frequency=frequency,
            significance_level=significance,
            min_depth=depth,
        )

    valid_haplotypes = findset(subconsensus_remapped, haplotype_tester)

    hapcalls = HaplotypeCall[]

    push!(hapcalls, HaplotypeCall(UInt64(0), 0.0, 0.0, 0.0, consensus_variant))

    for hap in valid_haplotypes
        length(hap) > 1 || continue
        push!(hapcalls, HaplotypeCall(hap, read_pool))
    end #for

    start_pos = min(leftposition.(fake_reads)...)
    end_pos = max(rightposition.(fake_reads)...)

    outdict = OrderedDict{String,Any}()
    outdict["version"] = HAPLINK_VERSION
    outdict["settings"] = args
    outdict["coverage"] = OrderedDict{String,Int}("start" => start_pos, "end" => end_pos)
    outdict["haplotypes"] = [_dict(h) for h in hapcalls]

    out_stream = isnothing(outfile) ? stdout : open(outfile, "w")
    YAML.write(out_stream, outdict)

    return 0
end #function

function _haplink_sequences(args::Dict{String,Any})
    reffile = args["reference"]
    hapfile = args["haplotypes"]
    outfile = args["output"]
    prefix = args["prefix"]

    # Read the haplotype file
    haplotype_input = YAML.load_file(hapfile)

    # Read the reference file
    refrecord = HapLink._first_record(reffile)
    refseq = FASTA.sequence(LongDNA{4}, refrecord)

    # Find the length of sequences
    startpos = haplotype_input["coverage"]["start"]
    endpos = haplotype_input["coverage"]["end"]

    # Parse the consensus variants into a Haplotype
    consensus_data = popfirst!(haplotype_input["haplotypes"])
    consensus_hap = Haplotype(refseq, consensus_data)
    consensus_seq = reconstruct(consensus_hap)
    consensus_ps = Pseudoread(startpos, endpos, consensus_hap)

    # Convert prefix to string
    prefix = isnothing(prefix) ? "" : prefix

    # Open the new fasta file
    FASTA.Writer(isnothing(outfile) ? stdout : open(outfile, "w")) do fasta_writer
        # Write the consensus sequence
        write(fasta_writer, FASTA.Record(consensus_ps; prefix=prefix, is_consensus=true))

        # Apply the variants of every other haplotype to a new haplotype
        for hap_dictionary in haplotype_input["haplotypes"]
            write(
                fasta_writer,
                FASTA.Record(
                    Pseudoread(
                        startpos, endpos, Haplotype(copy(consensus_seq), hap_dictionary)
                    );
                    prefix=prefix,
                    is_consensus=false,
                ),
            )
        end #for
    end #do

    return 0
end #function

@main
end #module
