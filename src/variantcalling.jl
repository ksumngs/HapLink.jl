"""
    phrederror(quality::Number)

Converts a PHRED33-scaled error number into the expected fractional error of basecall
"""
function phrederror(qual::Number)
    return 10^(-1 * qual / 10)
end #function

"""
    callvariants(snps::AbstractVector{SNP}, reads::AbstractPath, D_min::Int,
        Q_min::Int, x_min::Float64, f_min::Float64,
        α::Float64) where T <: Union{SAM.Record,BAM.Record}

Based on the aligned basecalls and stats in `bamcounts`, call variants.

# Arguments

  - `snps::AbstractVector{SNP}`: `Vector` of possible `SNPs` that variants will be called from
  - `reads::AbstractPath`: The path to an alignment file (SAM or BAM format) to call variants from
  - `D_min::Int`: minimum variant depth
  - `Q_min::Int`: minimum average PHRED-scaled quality at variant position
  - `x_min::Float64`: minimum average fractional distance from read end at variant position
  - `f_min::Float64`: minimum frequency of variant
  - `α::Float64`: significance level of variants by Fisher's Exact Test

# Returns

  - `Vector{VCF.Record}`: `VCF.Record`s that passed all the above filters
"""
function callvariants(
    snps::AbstractVector{<:SNP},
    reads::AbstractString,
    D_min::Int,
    Q_min::Int,
    x_min::Float64,
    f_min::Float64,
    α::Float64,
)

    # Setup a cache for the possible snps
    snpdata = DataFrame("snp" => snps)
    snpdata.location = ThreadsX.map(s -> location(s), snpdata.snp)

    # Calculate the depth
    snpdata.depth = ThreadsX.map(s -> depth(s, reads), snpdata.snp)

    # Remove snps that don't meet minimum depth
    filter!(:depth => d -> d >= D_min, snpdata)

    # Calculate the quality of remaining snps
    snpdata.quality = ThreadsX.map(s -> mean_quality(s, reads), snpdata.snp)

    # Remove snps that don't meet minimum quality
    filter!(:quality => q -> q >= Q_min, snpdata)

    # Position, converting [0 to 1] to [0 to 1 to 0]
    function pos_to_edge(x)
        if x <= 0.5
            return 2x
        else
            return 2 * (1 - x)
        end
    end #function

    # Calculate the position of remaining snps
    snpdata.position = ThreadsX.map(
        s -> pos_to_edge(mean_fractional_position(s, reads)), snpdata.snp
    )

    # Remove snps that are too close to the edge
    filter!(:position => p -> p >= x_min, snpdata)

    # Calculate the total depth at snp locations
    snpdata.totaldepth = ThreadsX.map(l -> depth(l, reads), snpdata.location)

    # Calculate the frequency of remaining snps
    snpdata.frequency = snpdata.depth ./ snpdata.totaldepth

    # Remove snps that aren't frequent enough
    filter!(:frequency => f -> f >= f_min, snpdata)

    # Calculate the depth of the reference base in
    snpdata.refdepth = ThreadsX.map(s -> depth(reference(s), reads), snpdata.snp)
    snpdata.totalqual = ThreadsX.map(l -> mean_quality(l, reads), snpdata.location)
    snpdata.ex_altdp = round.(Int, phrederror.(snpdata.totalqual) .* snpdata.totaldepth)
    snpdata.ex_refdp =
        round.(Int, (1 .- phrederror.(snpdata.totalqual)) .* snpdata.totaldepth)
    snpdata.pval =
        pvalue.(
            FisherExactTest.(
                snpdata.ex_altdp, snpdata.ex_refdp, snpdata.depth, snpdata.refdepth
            ),
        )

    # Calculate the Fisher's Exact probability of a variant base appearing due to sequencing
    # errors. This method is perfectly identical to the way iVar calls variants (See
    # https://github.com/andersen-lab/ivar/blob/v1.3.1/src/call_variants.cpp#L139), but
    # implemented in a slightly different way
    filter!(:pval => p -> p <= α, snpdata)

    # Return variant objects based on the remaining variant calls
    return VCF.Record.(eachrow(snpdata))
end #function

"""
    function savevcf(args...)

Save a VCF file populated with VCF records and metadata

# Arguments

  - `vcfs::AbstractVector{VCF.Record}`: `Vector` of `VCF.Record`s to write to file
  - `savepath::AbstractString`: path of the VCF file to write to. Will be overwritten
  - `refpath::AbstractString`: path of the reference genome used to call variants. The
    absolute path will be added to the `##reference` metadata
  - `D::Int`: mimimum variant depth used to filter variants. Will be added as `##FILTER`
    metadata
  - `Q::Number`: minimum PHRED quality used to filter variants. Will be added as `##FILTER`
    metadata
  - `x::Float64`: minimum fractional read position used to filter variants. Will be added as
    `##FILTER` metadata
  - `α::Float64`: Fisher's Exact Test significance level used to filter variants. Will be
    added as `##FILTER` metadata

Saves the variants in `vcfs` to a VCF file at `savepath`, adding the reference genome
`refpath`, the depth cutoff `D`, the quality cutoff `Q`, the position cutoff `x`, and the
significance cutoff `α` as metadata.
"""
function savevcf(
    vcfs::AbstractVector{VCF.Record},
    savepath::AbstractString,
    refpath::AbstractString,
    D::Int,
    Q::Number,
    x::Float64,
    α::Float64,
)

    # Convert read position to integer percent
    X = string(trunc(Int, x * 100))

    vheader = VCF.Header(
        VCF.MetaInfo.(
            split(
                """
                ##fileformat=VCFv4.2
                ##filedate=$(Dates.format(today(), "YYYYmmdd"))
                ##source=HapLink.jlv$(VERSION)
                ##reference=file://$(abspath(refpath))
                ##FILTER=<ID=d$D,Description="Variant depth below $D">
                ##FILTER=<ID=q$Q,Description="Quality below $Q">
                ##FILTER=<ID=x$X,Description="Position in outer $X% of reads">
                ##FILTER=<ID=sg,Description="Not significant at alpha=$α level by Fisher's Exact Test">
                ##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
                ##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate Depth">""",
                "\n",
            ),
        ),
        ["."],
    )

    f = VCF.Writer(open(savepath, "w"), vheader)

    for vcf in vcfs
        write(f, vcf)
    end #for

    return close(f)
end #function
