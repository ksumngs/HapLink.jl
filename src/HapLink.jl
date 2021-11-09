module HapLink

using FASTX
using HypothesisTests
using DataFrames
using FilePaths
using BioSequences

function main(args::Dict{String, Any})
    # 1. Analyze bam
    # 2. Call variants
    # 3. Call haplotypes

    bamfile        = args["bamfile"]
    reffile        = args["reference"]
    annotationfile = args["annotations"]
    Q_variant      = args["quality"]
    f_variant      = args["frequency"]
    x_variant      = args["position"]
    α_variant      = args["variant_significance"]
    D_variant      = args["variant_depth"]
    D_haplotype    = args["haplotype_depth"]

    bampath = Path(bamfile)
    prefix = isnothing(args["prefix"]) ? filename(bampath) : args["prefix"]

    variants = callvariants(analyzebam(bamfile, reffile),
        Q_variant, f_variant, x_variant, α_variant, D_variant)

    println(string(yamlize.(variants)...))
    println(string(vcfize.(variants)...))

end #function

function analyzebam(bamfile::String, reffile::String)
    bamanalysis = ""
    open(FASTA.Reader, reffile) do refreader
        for refrecord in refreader
            chromosome = FASTA.identifier(refrecord)
            seqlength = FASTA.seqlen(refrecord)
            bamanalysis = string(
                bamanalysis,
                readchomp(`bam-readcount -f $reffile $bamfile "$chromosome:1-$seqlength"`)
            )
        end #for
    end #do
    return string.(split(bamanalysis, "\n"))
end #function

function callvariants(bamcounts::Vector{String},
    minqual::Int, minfreq::Float64, minpos::Float64, α::Float64, D::Int)

    variantdata = bamcounts2bamdata(bamcounts)
    filter!(var -> var.base != var.reference_base, variantdata)
    filter!(var -> var.count >= D, variantdata)
    filter!(var -> var.avg_basequality >= minqual, variantdata)
    filter!(var -> var.avg_pos_as_fraction >= minpos, variantdata)
    filter!(var -> (var.count / var.depth) >= minfreq, variantdata)
    filter!(
        var -> pvalue(FisherExactTest(
            round(Int, phrederror(var.avg_basequality)*var.depth),
            round(Int, (1-phrederror(var.avg_basequality))*var.depth),
            var.count,
            var.depth
        )) <= α,
        variantdata
    )
    return Variant.(eachrow(variantdata))

end #function

function bamcounts2bamdata(bamcounts::AbstractVector{String})
    # Declare an empty bam stats data frame
    countsdata = DataFrame(
        chr                                  = String[],
        position                             = Int[],
        reference_base                       = String[],
        depth                                = Int[],
        base                                 = String[],
        count                                = Int[],
        avg_mapping_quality                  = Float64[],
        avg_basequality                      = Float64[],
        avg_se_mapping_quality               = Float64[],
        num_plus_strand                      = Int[],
        num_minus_strand                     = Int[],
        avg_pos_as_fraction                  = Float64[],
        avg_num_mismatches_as_fraction       = Float64[],
        avg_sum_mismatch_qualities           = Float64[],
        num_q2_containing_reads              = Int[],
        avg_distance_to_q2_start_in_q2_reads = Float64[],
        avg_clipped_length                   = Float64[],
        avg_distance_to_effective_3p_end     = Float64[]
    )

    # Transform the bam stats file
    for bamline in bamcounts
        # Split the base-independent stats by tabs
        bamfields = split(bamline, "\t")

        # Loop through the base-dependent stat blocks
        for i in 6:length(bamfields)
                # Split the base-dependent stats by colons
                basestats = split(bamfields[i], ":")

                # Parse the data into the correct types
                chr                                  = bamfields[1]
                position                             = parse(Int, bamfields[2])
                reference_base                       = bamfields[3]
                depth                                = parse(Int, bamfields[4])
                base                                 = basestats[1]
                count                                = parse(Int, basestats[2])
                avg_mapping_quality                  = parse(Float64, basestats[3])
                avg_basequality                      = parse(Float64, basestats[4])
                avg_se_mapping_quality               = parse(Float64, basestats[5])
                num_plus_strand                      = parse(Int, basestats[6])
                num_minus_strand                     = parse(Int, basestats[7])
                avg_pos_as_fraction                  = parse(Float64, basestats[8])
                avg_num_mismatches_as_fraction       = parse(Float64, basestats[9])
                avg_sum_mismatch_qualities           = parse(Float64, basestats[10])
                num_q2_containing_reads              = parse(Int, basestats[11])
                avg_distance_to_q2_start_in_q2_reads = parse(Float64, basestats[12])
                avg_clipped_length                   = parse(Float64, basestats[13])
                avg_distance_to_effective_3p_end     = parse(Float64, basestats[14])

                # Append the data to the dataframe
                push!(countsdata, [
                    chr,
                    position,
                    reference_base,
                    depth,
                    base,
                    count,
                    avg_mapping_quality,
                    avg_basequality,
                    avg_se_mapping_quality,
                    num_plus_strand,
                    num_minus_strand,
                    avg_pos_as_fraction,
                    avg_num_mismatches_as_fraction,
                    avg_sum_mismatch_qualities,
                    num_q2_containing_reads,
                    avg_distance_to_q2_start_in_q2_reads,
                    avg_clipped_length,
                    avg_distance_to_effective_3p_end
                ])
        end
    end

    return countsdata
end #function

function phrederror(qual::Number)
    return 10^(-1*qual/10)
end #function

struct Variant
    chromosome::String
    position::Int
    identifier::String
    referencebase::NucleotideSeq
    alternatebase::NucleotideSeq
    quality::Number
    filter::Symbol
    info::Dict{String,Any}
end #struct

function Variant(data::DataFrameRow)
    CHROM  = data.chr
    POS    = data.position
    ID     = "."
    QUAL   = data.avg_basequality
    FILTER = :PASS
    INFO   = Dict(
        "DP" => data.depth,
        "VD" => data.count
    )
    refbase    = data.reference_base
    altbase    = data.base

    # Check for insertion
    if first(altbase) == '+'
        altbase = string(refbase, altbase[2:end])
    end # if

    # Check for deletion
    if first(altbase) == '-'
        altbase = "-"
    end

    REF = LongDNASeq(refbase)
    ALT = LongDNASeq(altbase)

    return Variant(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
end #function

function yamlize(v::Variant)
    infostring = ""
    for n in v.info
        infostring = string(infostring, "      ", n[1], ": ", n[2], "\n")
    end
    return string(
        "  - chromosome: ",
        v.chromosome,
        "\n",
        "    position: ",
        string(v.position),
        "\n",
        "    identifier: ",
        string(v.identifier),
        "\n",
        "    referencebase: ",
        string(v.referencebase),
        "\n",
        "    alternatebase: ",
        string(v.alternatebase),
        "\n",
        "    quality: ",
        string(v.quality),
        "\n",
        "    filter: ",
        string(v.filter),
        "\n",
        "    info:\n",
        infostring
    )
end

function vcfize(v::Variant)
    return string(
        v.chromosome,             "\t",
        string(v.position),       "\t",
        v.identifier,             "\t",
        string(v.referencebase),  "\t",
        string(v.alternatebase),  "\t",
        string(trunc(v.quality)), "\t",
        string(v.filter),         "\t",
        join([string(n[1],"=",n[2]) for n in v.info], ";")
    )
end

end
