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
    _pos_to_edge(pos::Number)

Converts `pos` from a number between 0 (beginning) and 1 (end) to a number ranging from 0
(beginning) to 1 (middle) to 0 (end), i.e. convert a relative position to a relative
distance from the edge.
"""
function _pos_to_edge(pos::Number)
    pos > 0 && pos < 1 || error("`pos` must be between 0 and 1")
    return min(pos, 1 - pos)
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

"""
    _push_filter!(
        vc::VariationCall,
        label::Char,
        value::Union{Nothing,Number},
        filter::Function=(var, val) -> true,
    )

Adds a FILTER entry to `vc` of the form `"\$label\$value"` if `filter` returns true.

# Arguments
- `vc::VariationCall`: The [`VariationCall`](@ref) to annotate
- `label::Char`: The first character of the filter text to add
- `value::Union{Nothing,Number}`: The value to compare `vc` against, and to append to the
    filter text. If set to `nothing`, `_push_filter!` will return without evaluating or
    adding any filters, which may be useful for processing multiple inputs
- `filter::Function=(var, val) -> true`: A function handle to determine if a filter should
    be applied or not. Note that this function must return `true` **if and only if** the
    filter should be added to `vc`. The function will be passed `vc` and `value`. Defaults
    to always applying filters regardless of the values of `vc` and `value`.
"""
function _push_filter!(
    vc::VariationCall,
    label::Char,
    value::Union{Nothing,Number},
    filter::Function=(var, val) -> true,
)
    if isnothing(value)
        return nothing
    end #if

    filttext = "$label$value"
    if filter(vc, value)
        push!(vc.filter, filttext)
    end #if

    return filttext
end #function

"""
    call_variant(
        pileup::VariationPileup,
        α::Float64;
        D::Union{Nothing,Int}=nothing,
        Q::Union{Nothing,Float64}=nothing,
        X::Union{Nothing,Float64}=nothing,
        F::Union{Nothing,Float64}=nothing,
        S::Union{Nothing,Float64}=nothing,
    ) -> VariationCall

Calls variant from a `pileup`.

# Arguments
- `pileup::VariationPileup`: The pileup to call variants from
- `α::Float64`: Fisher's Exact Test significance (``α``) level to call variants

# Keywords

!!! note
    Leave any keyword undefined to skip filtering based on that field

- `D::Union{Nothing,Int}=nothing`: Minimum total read depth for a variant to be called
- `Q::Union{Nothing,Float64}=nothing`: Minimum average Phred quality for a variant to be
    called
- `X::Union{Nothing,Float64}=nothing`: Maximum average distance variant can be from read
    edge to be called
- `F::Union{Nothing,Float64}=nothing`: Minimum alternate frequency for a variant to be
    called
- `S::Union{Nothing,Float64}=nothing`: Maximum proportion of the number of times variant can
    appear on one strand versus the other
"""
function call_variant(
    pileup::VariationPileup,
    α::Float64;
    D::Union{Nothing,Int}=nothing,
    Q::Union{Nothing,Float64}=nothing,
    X::Union{Nothing,Float64}=nothing,
    F::Union{Nothing,Float64}=nothing,
    S::Union{Nothing,Float64}=nothing,
)
    vc = VariationCall(pileup)
    _push_filter!(vc, 'd', D, (var, val) -> altdepth(var) < val)
    _push_filter!(vc, 'q', Q, (var, val) -> quality(var) < val)
    _push_filter!(vc, 'x', X, (var, val) -> _pos_to_edge(readpos(var)) < val)
    _push_filter!(vc, 'f', F, (var, val) -> frequency(var) < val)
    _push_filter!(vc, 's', S, (var, val) -> _pos_to_edge(strand_bias(var)) < val)
    _push_filter!(vc, 'a', α, (var, val) -> p_value(var) > val)

    if isempty(filters(vc))
        push!(vc.filter, "PASS")
    end #if

    return vc
end #function
