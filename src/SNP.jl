using BioSymbols
using GenomicFeatures

struct SNP{S, T <: NucleicAcid}
    location::Interval{S}
    refbase::T
    altbase::T
end #struct
