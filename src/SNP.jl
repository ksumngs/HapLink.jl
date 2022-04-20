using BioGenerics
using BioSymbols
using GenomicFeatures
using Lazy
using VariantCallFormat
using XAM

export SNP
export altbase
export location
export refbase
export reference

struct SNP{S<:NucleicAcid}
    location::Interval
    refbase::S
    altbase::S
end #struct

"""
    SNP(chrom::String, pos::Int64, ref::S, alt::S) where S<:NucleicAcid

A single-nuclotide polymorphism (SNP) located on the chromosome `chrom` at 1-based sequence
postion `pos`, that alters `ref` into `alt`.
"""
function SNP(chrom::String, pos::Int64, ref::S, alt::S) where {S<:NucleicAcid}
    return SNP(Interval(chrom, pos, pos), ref, alt)
end #function

"""
    SNP(v::VCF.Record)

Create an object to represent the mutation in `v` without metadata
"""
function SNP(v::VCF.Record)
    # Check if the variant is a multi-insertion
    if length(VCF.alt(v)) > 1 || length(first(VCF.alt(v))) > 1
        throw(
            ArgumentError(
                "VCF record cannot be converted to a SNP because it has more than one alternate base",
            ),
        )
    end #if

    # Check if the bases are RNA
    nuctype = DNA
    if 'U' in VCF.ref(v) || 'U' in first(VCF.alt(v))
        nuctype = RNA
    end

    # Create the new SNP object
    return SNP(
        Interval(VCF.chrom(v), VCF.pos(v), VCF.pos(v)),
        nuctype(first(VCF.ref(v))),
        nuctype(first(first(VCF.alt(v)))),
    )
end #struct

function location(s::SNP)
    return s.location
end #function

function refbase(s::SNP)
    return s.refbase
end #function

function altbase(s::SNP)
    return s.altbase
end #function

# Forward interval methods
@forward SNP.location BioGenerics.seqname,
BioGenerics.leftposition, BioGenerics.rightposition,
BioGenerics.metadata

function Base.isless(a::SNP, b::SNP)
    return location(a) < location(b) && refbase(a) < refbase(b) && altbase(a) < altbase(b)
end #function

function Base.show(io::IO, s::SNP)
    return print(io, "$(seqname(s)):$(leftposition(s)) ($(refbase(s)) -> $(altbase(s)))")
end #function

"""
    reference(s::SNP)

Create a `SNP` representing the reference base of `s`.
"""
function reference(s::SNP)
    return SNP(s.location, s.refbase, s.refbase)
end #function

function VariantCallFormat.VCF.Record(
    s::SNP;
    id::Union{String,Nothing}=nothing,
    qual::Union{Float64,Nothing}=nothing,
    filter::String="PASS",
    depth::Union{Int64,Nothing}=nothing,
    altdepth::Union{Int64,Nothing}=nothing,
    info::Dict{<:Any,<:Any}=Dict(),
)

    # Convert the id and quality
    id_str = isnothing(id) ? "." : copy(id)
    qual_str = isnothing(qual) ? "." : trunc(qual; digits=1)

    # Copy the dictionary over
    info_dict = copy(info)

    # Add depth and altdepth to the dictionary
    if !isnothing(depth)
        info_dict["DP"] = depth
    end #if
    if !isnothing(altdepth)
        info_dict["AD"] = altdepth
    end #if

    # Convert the snp and metadata to string, then to VCF.Record
    return VCF.Record(
        join(
            [
                seqname(s),
                leftposition(s),
                id_str,
                refbase(s),
                altbase(s),
                qual_str,
                filter,
                join(["$(n[1])=$(n[2])" for n in info_dict], ";"),
            ],
            "\t",
        ),
    )
end #function

function VariantCallFormat.VCF.Record(
    snp::SNP,
    reads::AbstractVector{T};
    id::Union{String,Nothing}=nothing,
    filter::String="PASS",
) where {T<:Union{SAM.Record,BAM.Record}}
    return VCF.Record(
        snp;
        id=id,
        qual=mean_quality(snp, reads),
        filter=filter,
        depth=depth(location(snp), reads),
        altdepth=depth(snp, reads),
    )
end #function

function VariantCallFormat.VCF.Record(record::DataFrameRow)
    return VCF.Record(
        record.snp; qual=record.quality, altdepth=record.depth, depth=record.totaldepth
    )
end #function
