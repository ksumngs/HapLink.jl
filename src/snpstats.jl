using XAM

export doescontain

"""
    doescontain(snp::SNP, rec::Union{SAM.Record,BAM.Record})

Determines if the location and base mutation of `snp` are present in the sequence of `rec`

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> # SAM spec example record, with an added mutation at postion 12

julia> samrecord = SAM.Record("r001\\t99\\tref\\t7\\t30\\t8M2I4M1D3M\\t=\\t37\\t39\\tTTAGACAAAGGATACTG\\t*");

julia> doescontain(SNP("ref", 6, DNA_T, DNA_C), samrecord)
true

julia> doescontain(SNP("ref", 17, DNA_T, DNA_A), samrecord)
false
```
"""
function doescontain(snp::SNP, rec::Union{SAM.Record,BAM.Record})
    return doescontain(snp.location, rec) &&
           first(basesat(snp.location, rec)) == snp.altbase
end #function
