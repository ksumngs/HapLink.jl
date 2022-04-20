export countbasestats
export possible_snps
export transformbamcounts

"""
    possible_snps(reference::AbstractVector{FASTA.Record})

Get an array of every possible `SNP` that could be mutated from `reference`.

# Example

```jldoctest
julia> using BioSymbols, FASTX, GenomicFeatures

julia> first(possible_snps([FASTA.Record(">ref\\nAGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT")]), 2)
2-element Vector{SNP}:
 ref:1 (A -> C)
 ref:1 (A -> G)
```
"""
function possible_snps(reference::AbstractVector{FASTA.Record})
    snps = SNP[]
    for chrom in reference
        for i in 1:FASTA.seqlen(chrom)
            alph = typeof(FASTA.sequence(chrom)[i]) == DNA ? ACGT : ACGU
            for N in alph
                if N != FASTA.sequence(chrom)[i]
                    push!(
                        snps, SNP(FASTA.identifier(chrom), i, FASTA.sequence(chrom)[i], N)
                    )
                end #if
            end #for
        end #for
    end #for
    return snps
end #function
