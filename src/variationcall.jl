"""
    VariationCall

Represents a `Variation` that is supported by aligned reads with sufficient metadata to
support or refute its validity. It is designed to be converted into a line in [Variant Call
Format](https://github.com/samtools/hts-specs#variant-calling-data-files), or a
[`VCF.Record`](https://rasmushenningsson.github.io/VariantCallFormat.jl/stable/vcf-bcf/).

# Fields
- `variation::Variation`: The `Variation` of this call
- `quality::Union{Nothing,Number}`: Phred quality of the basecall for `variation`
- `filter::Union{Nothing,Vector{String}}`: Indicator if `variation` has passed filters and
  is actually a variant call, or a list of criteria that have failed it
- `depth::Union{Nothing,UInt}`: The number of reads that cover `leftposition(variation)`
- `strandbias::Union{Nothing,Float64}`: The fraction of times `variation` appears on a
  positive strand
- `altdepth::Union{Nothing,UInt}`: The number of types `variation` occurs
- `readpos::Union{Nothing,UInt}`: The averagerelative position of `variation` in each read
- `pvalue::Union{Nothing,Float64}`: The Fisher's Exact Test ``p``-value of this call
"""
Base.@kwdef struct VariationCall
    variation::Variation
    quality::Union{Nothing,Number} = nothing
    filter::Union{Nothing,Vector{String}} = nothing
    depth::Union{Nothing,UInt} = nothing
    strandbias::Union{Nothing,Float64} = nothing
    altdepth::Union{Nothing,UInt} = nothing
    readpos::Union{Nothing,Float64} = nothing
    pvalue::Union{Nothing,Float64} = nothing
end #struct

function VariationCall(vp::VariationPileup)
    v = variation(vp)
    q = mean(quality(vp))
    f = String[]
    d = depth(vp)
    s = strand_bias(vp)
    a = UInt(altdepth(vp))
    x = mean(readpos(vp))
    p = variation_test(Int(d), Int(a), q)

    return VariationCall(v, q, f, d, s, a, x, p)
end #function

variation(vc::VariationCall) = vc.variation
quality(vc::VariationCall) = vc.quality
filters(vc::VariationCall) = vc.filter
depth(vc::VariationCall) = vc.depth
strand_bias(vc::VariationCall) = vc.strandbias
altdepth(vc::VariationCall) = vc.altdepth
readpos(vc::VariationCall) = vc.readpos
p_value(vc::VariationCall) = vc.pvalue

frequency(vc::VariationCall) = altdepth(vc) / depth(vc)

"""
    _phrederror(quality::Number)

Converts a PHRED33-scaled error number into the expected fractional error of basecall
"""
function _phrederror(qual::Number)
    return 10^(-1 * qual / 10)
end #function

"""
    variation_test(depth::Int, altdepth::Int, quality::Float64)

Conducts a Fisher's Exact Test to deterimine the likelihood of a variant with total `depth`
and variation depth `altdepth` occuring, given an average basecall `quality`. Returns the
``p``-value of the test.
"""
function variation_test(depth::Int, altdepth::Int, quality::Float64)
    refdepth = depth - altdepth
    expected_altdepth = round(Int, _phrederror(quality) * depth)
    expected_refdepth = round(Int, (1 - _phrederror(quality)) * depth)
    return pvalue(FisherExactTest(expected_altdepth, expected_refdepth, depth, refdepth))
end #function
