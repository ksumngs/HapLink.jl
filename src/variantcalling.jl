export callvariants

"""
    phrederror(quality::Number)

Converts a PHRED33-scaled error number into the expected fractional error of basecall
"""
function phrederror(qual::Number)
    return 10^(-1 * qual / 10)
end #function

"""
    callvariants(snps::AbstractVector{SNP}, reads::AbstractVector{T}, D_min::Int,
        Q_min::Int, x_min::Float64, f_min::Float64,
        α::Float64) where T <: Union{SAM.Record,BAM.Record}

Based on the aligned basecalls and stats in `bamcounts`, call variants.

# Arguments
- `snps::AbstractVector{SNP}`: `Vector` of possible `SNPs` that variants will be called from
- `reads::AbstractVector{T} where T <: Union{SAM.Record,BAM.Record}`: A set of aligned
    sequencing reads from which to verify variant calls
- `D_min::Int`: minimum variant depth
- `Q_min::Int`: minimum average PHRED-scaled quality at variant position
- `x_min::Float64`: minimum average fractional distance from read end at variant position
- `f_min::Float64`: minimum frequency of variant
- `α::Float64`: significance level of variants by Fisher's Exact Test

# Returns
- `Vector{SNP}`: [`SNP`](@ref)s that passed all the above filters
"""
function callvariants(
    snps::AbstractVector{<:SNP},
    reads::AbstractVector{T},
    D_min::Int,
    Q_min::Int,
    x_min::Float64,
    f_min::Float64,
    α::Float64,
) where {T <: Union{SAM.Record, BAM.Record}}

    # Ensure we don't mutate the input reads
    snpdata = copy(snps)

    # Simple metric filters
    # Depth
    filter!(s -> depth(s, reads) >= D_min, snpdata)

    # Quality
    filter!(s -> mean_quality(s, reads) >= Q_min, snpdata)

    # Position, converting 0 to 1 to 0 to 1 to 0
    function pos_to_edge(x)
        if x <= 0.5
            return 2x
        else
            return 2*(1-x)
        end
    end #function
    filter!(s -> pos_to_edge(fractional_position(s, reads)) >= x_min, snpdata)

    filter!(s -> frequency(s, reads) >= f_min, snpdata)

    # Calculate the Fisher's Exact probability of a variant base appearing due to sequencing
    # errors. This method is perfectly identical to the way iVar calls variants (See
    # https://github.com/andersen-lab/ivar/blob/v1.3.1/src/call_variants.cpp#L139), but
    # implemented in a slightly different way
    filter!(
        s ->
            pvalue(
                FisherExactTest(
                    round(Int, phrederror(mean_quality(s, reads)) * depth(s.location, reads)),
                    round(Int, (1 - phrederror(mean_quality(s, reads))) * depth(s.location, reads)),
                    depth(s, reads),
                    depth(reference(s), reads),
                ),
            ) <= α,
        snpdata,
    )

    # Return variant objects based on the remaining variant calls
    return snpdata
end #function

"""
    savevcf(vars::AbstractVector{Variant}, savepath::String, refpath::String, D::Int,
        Q::Number, x::Float64, α::Float64)

Save a VCF file populated with `vars`

# Arguments
- `vars::AbstractVector{Variant}`: `Vector` of [`Variant`](@ref)s to write to file
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

Saves the variants in `vars` to a VCF file at `savepath`, adding the reference genome
`refpath`, the depth cutoff `D`, the quality cutoff `Q`, the position cutoff `x`, and the
significance cutoff `α` as metadata.
"""
function savevcf(
    vars::AbstractVector{Variant},
    savepath::AbstractString,
    refpath::AbstractString,
    D::Int,
    Q::Number,
    x::Float64,
    α::Float64,
)

    # Convert read position to integer percent
    X = string(trunc(Int, x * 100))

    # Open the file via clobbering
    open(savepath, "w") do f
        # Write headers
        write(f, "##fileformat=VCFv4.2\n")
        write(f, string("##filedate=", Dates.format(today(), "YYYYmmdd"), "\n"))
        write(f, string("##source=HapLink.jlv", VERSION, "\n"))
        write(f, string("##reference=file://", abspath(refpath), "\n"))

        # Write filter metadata
        write(f, "##FILTER=<ID=d$D,Description=\"Variant depth below $D\">\n")
        write(f, "##FILTER=<ID=q$Q,Description=\"Quality below $Q\">\n")
        write(f, "##FILTER=<ID=x$X,Description=\"Position in outer $X% of reads\">\n")
        write(
            f,
            "##FILTER=<ID=sg,Description=\"Not significant at α=$α level by Fisher's Exact Test\">\n",
        )

        # Add descriptions of the info tags I chose to include
        # TODO: Find a way for these _not_ to be hard-coded in here
        write(f, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        write(f, "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Alternate Depth\">\n")

        # Write the header line?
        # TODO: Check this header against VCF spec
        write(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Write every variant out
        for var in vars
            write(f, string(serialize_vcf(var), "\n"))
        end #for
    end #do
end #function
