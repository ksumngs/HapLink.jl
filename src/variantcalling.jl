export callvariants

"""
    phrederror(quality::Number)

Converts a PHRED33-scaled error number into the expected fractional error of basecall
"""
function phrederror(qual::Number)
    return 10^(-1 * qual / 10)
end #function

"""
    callvariants(bamcounts::AbstractDataFrame, D_min::Int, Q_min::Int, x_min::Float64,
        f_min::Float64, α::Float64)

Based on the aligned basecalls and stats in `bamcounts`, call variants.

# Arguments
- `bamcounts::AbstractDataFrame`: `DataFrame` containing the output from `bam-readcount`
- `D_min::Int`: minimum variant depth
- `Q_min::Int`: minimum average PHRED-scaled quality at variant position
- `x_min::Float64`: minimum average fractional distance from read end at variant position
- `f_min::Float64`: minimum frequency of variant
- `α::Float64`: significance level of variants by Fisher's Exact Test

# Returns
- `Vector{Variant}`: [`Variant`](@ref)s that passed all the above filters
"""
function callvariants(
    bamcounts::AbstractDataFrame,
    D_min::Int,
    Q_min::Int,
    x_min::Float64,
    f_min::Float64,
    α::Float64,
)

    # Ensure we don't mutate the input reads
    variantdata = copy(bamcounts)

    # Exclude any non-variants i.e. reference bases
    filter!(var -> var.base != var.reference_base, variantdata)

    # Simple metric filters
    filter!(var -> var.count >= D_min, variantdata)
    filter!(var -> var.avg_basequality >= Q_min, variantdata)
    filter!(var -> var.avg_pos_as_fraction >= x_min, variantdata)
    filter!(var -> (var.count / var.depth) >= f_min, variantdata)

    # Calculate the Fisher's Exact probability of a variant base appearing due to sequencing
    # errors. This method is perfectly identical to the way iVar calls variants (See
    # https://github.com/andersen-lab/ivar/blob/v1.3.1/src/call_variants.cpp#L139), but
    # implemented in a slightly different way
    filter!(
        var ->
            pvalue(
                FisherExactTest(
                    round(Int, phrederror(var.avg_basequality) * var.depth),
                    round(Int, (1 - phrederror(var.avg_basequality)) * var.depth),
                    var.count,
                    var.depth,
                ),
            ) <= α,
        variantdata,
    )

    # Return variant objects based on the remaining variant calls
    return Variant.(eachrow(variantdata))
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
