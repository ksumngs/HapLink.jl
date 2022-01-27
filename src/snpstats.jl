using XAM

export doescontain
export depth

"""
    doescontain(snp::SNP, rec::Union{SAM.Record,BAM.Record})

Determines if the location and base mutation of `snp` are present in the sequence of `rec`

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> # SAM spec example record, with an added mutation at postion 12

julia> samrecord = SAM.Record("r001\\t99\\tref\\t7\\t30\\t8M2I4M1D3M\\t=\\t37\\t39\\tTTAGACAAAGGATACTG\\t*");

julia> doescontain(SNP("ref", 12, DNA_T, DNA_C), samrecord)
true

julia> doescontain(SNP("ref", 17, DNA_T, DNA_A), samrecord)
false
```
"""
function doescontain(snp::SNP, rec::Union{SAM.Record,BAM.Record})
    return doescontain(snp.location, rec) &&
           first(basesat(snp.location, rec)) == snp.altbase
end #function

"""
    depth(snp::SNP, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Calculate the number of times `snp` is represented within `reads`

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> # All SAM spec example sequences, with an added mutation at position 17 in reads 1 and 4

julia> samrecords = SAM.Record.([
           "r001\\t99\\tref\\t7\\t30\\t8M2I4M1D3M\\t=\\t37\\t39\\tTTAGATAAAGGACACTG\\t*",
           "r002\\t0\\tref\\t9\\t30\\t3S6M1P1I4M\\t*\\t0\\t0\\tAAAAGATAAGGATA\\t*",
           "r003\\t0\\tref\\t9\\t30\\t5S6M\\t*\\t0\\t0\\tGCCTAAGCTAA\\t*\\tSA:Z:ref,29,-,6H5M,17,0;",
           "r004\\t0\\tref\\t16\\t30\\t6M14N5M\\t*\\t0\\t0\\tACAGCTTCAGC\\t*",
           "r003\\t2064\\tref\\t29\\t17\\t6H5M\\t*\\t0\\t0\\tTAGGC\\t*\\tSA:Z:ref,9,+,5S6M,30,1;",
           "r001\\t147\\tref\\t37\\t30\\t9M\\t=\\t7\\t-39\\tCAGCGGCAT\\t*\\tNM:i:1"
       ]);

julia> depth(SNP("ref", 17, DNA_T, DNA_C), samrecords)
2
```
"""
function depth(snp::SNP, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}
    return count(r -> doescontain(snp, r), reads)
end #function
