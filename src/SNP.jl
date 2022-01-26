using BioSymbols
using GenomicFeatures
using VariantCallFormat

export SNP

struct SNP{S<:NucleicAcid}
    location::Interval
    refbase::S
    altbase::S
end #struct

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
