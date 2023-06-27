"""
    VariationPileup

Summarizes the basecalls that support a `Variation` within a set of alignments.

# Fields
- `variation::Variation`: The Variation of this pileup entry
- `depth::UInt`: The number of times the position of `variation` appears in this set of
  alignments
- `readpos::Vector{Float64}`: The relative positions of `variation` within each read
- `quality::Vector{Float64}`: The phred quality of `variation` within each read
- `strand::Vector{Strand}`: Which strand each `variation` is found on
"""
Base.@kwdef struct VariationPileup
    variation::Variation
    depth::UInt
    readpos::Vector{Float64} = Float64[]
    quality::Vector{Float64} = Float64[]
    strand::Vector{Strand} = Strand[]
end #struct

"""
    variation(vp::VariationPileup)

Gets the `Variation` of `vp`
"""
variation(vp::VariationPileup) = vp.variation

"""
    depth(vp::VariationPileup)

Gets the number of times the position of `vp.variation` appears total (variant and
wild-type)
"""
depth(vp::VariationPileup) = vp.depth

"""
    readpos(vp::VariationPileup)

Gets the relative positions of `vp.variation`
"""
readpos(vp::VariationPileup) = vp.readpos

"""
    quality(vp::VariationPileup)

Gets the phred qualities of `vp.variation`
"""
quality(vp::VariationPileup) = vp.quality

"""
    strand(vp::VariationPileup)

Gets the strands of `vp.variation`
"""
strand(vp::VariationPileup) = vp.strand

"""
    altdepth(vp::VariationPileup)

Gets the number of times `vp.variation` appears
"""
function altdepth(vp::VariationPileup)
    return min(length(readpos(vp)), length(quality(vp)), length(strand(vp)))
end #function

"""
    frequency(vp::VariationPileup)

Gets the alternate allele frequency of `vp.variation`
"""
frequency(vp::VariationPileup) = altdepth(vp) / depth(vp)

"""
    strand_bias(vp::VariationPileup)

Gets the frequency of positive strands that `variation` appears on relative to all
`variation` reads
"""
function strand_bias(vp::VariationPileup)
    return count(s -> s == STRAND_POS, strand(vp)) / altdepth(vp)
end #function

function Base.:(==)(x::VariationPileup, y::VariationPileup)
    return variation(x) == variation(y) &&
           depth(x) == depth(y) &&
           readpos(x) == readpos(y) &&
           quality(x) == quality(y) &&
           strand(x) == strand(y)
end #function

function Base.hash(x::VariationPileup, h::UInt)
    return hash(
        VariationPileup,
        hash((variation(x), depth(x), readpos(x), quality(x), strand(x)), h),
    )
end

function _depth_dict(sam::Union{AbstractString,AbstractPath}, ref::FASTA.Record)
    dd = Dict{Interval,UInt}()

    rname = FASTA.identifier(ref)
    rlen = length(FASTA.sequence(ref))
    sams = last(_indexed_reader(sam, Interval(rname, 1, rlen)))

    for i in 1:rlen
        int = Interval(rname, i, i)
        dd[int] = 0
    end #for

    XAM_LEFTPOSITION = _issam(sam) ? SAM.position : BAM.position
    XAM_RIGHTPOSITION = _issam(sam) ? SAM.rightposition : BAM.rightposition

    for s in sams
        pos = XAM_LEFTPOSITION(s):XAM_RIGHTPOSITION(s)

        for p in pos
            int = Interval(rname, p, p)
            dd[int] += 1
        end #for
    end #for

    return dd
end #function

"""
    pileup(sam::AbstractPath, ref::AbstractPath) -> Vector{VariationPileup}

Generates a pileup of `Variation`s based on the alignments in `sam` aligned against `ref`.
"""
function pileup(
    sam::Union{AbstractString,AbstractPath}, ref::Union{AbstractString,AbstractPath}
)
    vis = variationinfos(sam, ref)
    vip = Dict{Variation,VariationPileup}()

    rrec = _first_record(ref)
    rname = FASTA.identifier(rrec)

    dd = _depth_dict(sam, rrec)

    vs = Variation[]
    for vi in vis
        push!(vs, variation(vi))
    end #for
    unique!(vs)

    for v in vs
        int = Interval(rname, leftposition(v), leftposition(v))
        d = dd[int]
        vip[v] = VariationPileup(; variation=v, depth=d)
    end #for

    for vi in vis
        v = variation(vi)
        push!(vip[v].readpos, readpos(vi))
        push!(vip[v].quality, quality(vi))
        push!(vip[v].strand, strand(vi))
    end #for

    return collect(values(vip))
end #function
