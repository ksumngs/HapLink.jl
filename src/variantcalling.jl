"""
    phrederror(quality::Number)

Converts a PHRED33-scaled error number into the expected fractional error of basecall
"""
function phrederror(qual::Number)
    return 10^(-1*qual/10)
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
    α::Float64
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
        var -> pvalue(FisherExactTest(
            round(Int, phrederror(var.avg_basequality)*var.depth),
            round(Int, (1-phrederror(var.avg_basequality))*var.depth),
            var.count,
            var.depth
        )) <= α,
        variantdata
    )

    # Return variant objects based on the remaining variant calls
    return Variant.(eachrow(variantdata))

end #function
