using BioGenerics
using BioSymbols
using GenomicFeatures
using Lazy
using VariantCallFormat

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
function SNP(chrom::String, pos::Int64, ref::S, alt::S) where S<:NucleicAcid
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
@forward SNP.location BioGenerics.seqname, BioGenerics.leftposition, BioGenerics.rightposition, BioGenerics.metadata

function Base.isless(a::SNP, b::SNP)
    return location(a) < location(b) && refbase(a) < refbase(b) && altbase(a) < altbase(b)
end #function

function Base.show(io::IO, s::SNP)
    print(io, "$(seqname(s)):$(leftposition(s)) ($(refbase(s)) -> $(altbase(s)))")
end #function

"""
    reference(s::SNP)

Create a `SNP` representing the reference base of `s`.
"""
function reference(s::SNP)
    return SNP(s.location, s.refbase, s.refbase)
end #function
