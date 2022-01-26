using BioSymbols
using GenomicFeatures

struct SNP{S<:NucleicAcid}
    location::Interval
    refbase::S
    altbase::S
end #struct
