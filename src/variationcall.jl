"""
    VariationCall

Represents a `Variation` that is supported by aligned reads with sufficient metadata to
support or refute its validity. It is designed to be converted into a line in [Variant Call
Format](https://github.com/samtools/hts-specs#variant-calling-data-files), or a
[`VCF.Record`](https://rasmushenningsson.github.io/VariantCallFormat.jl/stable/vcf-bcf/).

# Fields
- `variation::Variation`: The `Variation` of this call
- `quality::Union{Nothing,Number}`: Phred quality of the basecall for `variation`
- `filter::Vector{String}`: Indicator if `variation` has passed filters and is actually
  a variant call, or a list of criteria that have failed it
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
    filter::Vector{String} = String[]
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

"""
    variation(vc::VariationCall) -> SequenceVariation.Variation

Gets the `Variation` of a [`VariationCall`](@ref)
"""
variation(vc::VariationCall) = vc.variation

"""
    quality(vc::VariationCall) -> Union{Nothing,Float64}

Gets the average phred quality score of `vc`, if known. Returns `nothing` if unknown.

See also [`quality(::VariationInfo)`](@ref), [`quality(::VariationPileup)`](@ref)
"""
quality(vc::VariationCall) = vc.quality

"""
    filters(vc::VariationCall) -> Vector{String}

Gets all filters that have been applied to `vc`. Note that an empty FILTER entry is not
permitted under the [VCF
spec](https://github.com/samtools/hts-specs#variant-calling-data-files), and an empty array
should not automatically be considered to have `PASS`ed all filters.
"""
filters(vc::VariationCall) = vc.filter

"""
    depth(vc::VariationCall) -> Union{Nothing,UInt}

Gets the number of times the position of `vc` appears total. Returns `nothing` if unknown.

See also [`depth(::VariationPileup)`](@ref)
"""
depth(vc::VariationCall) = vc.depth

"""
    strand_bias(vc::VariationCall) -> Union{Nothing,Float64}

Gets the fraction of times `vc` appears on the positive strand. Returns `nothing` if
unknown.

See also [`strand(::VariationInfo)`](@ref), [`strand(::VariationPileup)`](@ref)
"""
strand_bias(vc::VariationCall) = vc.strandbias

"""
    altdepth(vc::VariationCall) -> Union{Nothing,UInt}

Gets the number of times `vc` appears. Returns `nothing` if unknown.

See also [`altdepth(::VariationPileup)`](@ref)
"""
altdepth(vc::VariationCall) = vc.altdepth

"""
    readpos(vc::VariationCall) -> Union{Nothing,Float64}

Gets the average relative position of `vc`. Returns `nothing` if unknown.

See also [`readpos(::VariationPileup)`](@ref)
"""
readpos(vc::VariationCall) = vc.readpos

"""
    p_value(vc::VariationCall) -> Union{Nothing,Float64}

Gets the ``p``-value of the observed statistic of `vc`. Returns `nothing` if unknown.

See also [`variation_test`](@ref)
"""
p_value(vc::VariationCall) = vc.pvalue

"""
    frequency(vc::VariationCall) -> Float64

Gets the alternate allele frequency of `vc`

See also [`frequency(::VariationPileup)`](@ref)
"""
frequency(vc::VariationCall) = altdepth(vc) / depth(vc)

"""
    vcf(vc::VariationCall, refname::AbstractString) -> VCF.Record

Converts `vc` into a `VCF.Record`. `refname` is required and used as the `CHROM` field in
the record.
"""
function vcf(vc::VariationCall, refname::AbstractString)
    chrom = refname
    pos = leftposition(variation(vc))
    id = "."
    ref = refbases(variation(vc))
    alt = altbases(variation(vc))

    filt = join(filters(vc), ";")

    qual = isnothing(quality(vc)) ? "." : quality(vc)

    info = Dict()
    info["DP"] = isnothing(depth(vc)) ? "." : depth(vc)
    info["SB"] = isnothing(strand_bias(vc)) ? "." : strand_bias(vc)
    info["AD"] = isnothing(altdepth(vc)) ? "." : altdepth(vc)
    info["AF"] = isnothing(frequency(vc)) ? "." : frequency(vc)
    info["RP"] = isnothing(readpos(vc)) ? "." : readpos(vc)
    info["PVAL"] = isnothing(p_value(vc)) ? "." : p_value(vc)

    info_string = join(["$k=$v" for (k, v) in info], ";")

    vcf_string = join([chrom, pos, id, ref, alt, qual, filt, info_string], "\t")

    return VCF.Record(vcf_string)
end #function

function _vcf_header(
    refpath::Union{AbstractPath,AbstractString},
    α::Float64;
    D::Union{Nothing,Int}=nothing,
    Q::Union{Nothing,Float64}=nothing,
    X::Union{Nothing,Float64}=nothing,
    F::Union{Nothing,Float64}=nothing,
    S::Union{Nothing,Float64}=nothing,
)
    MI = VCF.MetaInfo

    abs_ref = absolute(Path(refpath))

    _to_percent(x) = isnothing(x) ? nothing : string(trunc(Int, x * 100))
    x = _to_percent(X)
    f = _to_percent(F)
    s = _to_percent(S)

    metas = MI[]

    push!(metas, MI("##fileformat=VCFv4.3"))
    push!(metas, MI("##filedate=$(Dates.format(today(), "YYYmmdd"))"))
    push!(metas, MI("##source=HapLink.jlv$VERSION"))
    push!(metas, MI("##reference=file://$abs_ref"))
    push!(
        metas,
        MI(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Combined depth across samples\">",
        ),
    )
    push!(metas, MI("##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias\">"))
    push!(
        metas,
        MI(
            "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Total read depth for each allele\">",
        ),
    )
    push!(
        metas,
        MI(
            "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency for each ALT allele\">",
        ),
    )
    push!(
        metas,
        MI("##INFO=<ID=RP,Number=1,Type=Float,Description=\"Read position of allele\">"),
    )
    push!(
        metas,
        MI(
            "##INFO=<ID=PVAL,Number=1,Type=Float,Description=\"Fisher's exact p-value of variant call\">",
        ),
    )
    push!(
        metas,
        MI(
            "##FILTER=<ID=a$α,Description=\"Not significant at alpha=$α level by Fisher's Exact Test\">",
        ),
    )
    isnothing(D) ||
        push!(metas, MI("##FILTER=<ID=d$D,Description=\"Variant depth below $D\">"))
    isnothing(Q) || push!(metas, MI("##FILTER=<ID=q$Q,Description=\"Quality below $Q\">"))
    isnothing(x) ||
        push!(metas, MI("##FILTER=<ID=x$x,Description=\"Position in outer $x% of reads\">"))
    isnothing(f) ||
        push!(metas, MI("##FILTER=<ID=f$f,Description=\"Frequency below $f%\">"))
    isnothing(s) || push!(metas, "##FILTER=<ID=s$s,Description=\"Strand bias above $s%\">")

    return VCF.Header(metas, ["."])
end #function

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
    pos >= 0 && pos <= 1 || error("`pos` must be between 0 and 1")
    return pos <= 0.5 ? pos * 2 : (1 - pos) * 2
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
