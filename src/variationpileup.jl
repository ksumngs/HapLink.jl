Base.@kwdef struct VariationPileup
    variation::Variation
    depth::UInt
    readpos::Vector{Float64} = Float64[]
    quality::Vector{Float64} = Float64[]
    strand::Vector{Strand} = Strand[]
end #struct

variation(vp::VariationPileup) = vp.variation
depth(vp::VariationPileup) = vp.depth
readpos(vp::VariationPileup) = vp.readpos
quality(vp::VariationPileup) = vp.quality
strand(vp::VariationPileup) = vp.strand
function altdepth(vp::VariationPileup)
    return min(length(readpos(vp)), length(quality(vp)), length(strand(vp)))
end #function
frequency(vp::VariationPileup) = altdepth(vp) / depth(vp)

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
    rlen = FASTA.seqlen(ref)
    sams = last(_indexed_reader(sam, Interval(rname, 1, rlen)))

    for i in 1:FASTA.seqlen(ref)
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

    return values(vip)
end #function
