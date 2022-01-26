using BioSymbols
using GenomicFeatures

struct SNP{S <: NucleicAcid}
    location::Interval
    refbase::T
    altbase::T
end #struct
