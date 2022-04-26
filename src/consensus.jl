"""
    consensus(refseq::LongDNASeq, vars::AbstractArray{Variant}; freq::Float64=0.5)

Generate the consensus sequence of `refseq` when mutated by `vars`, excluding any items from
`vars` that have a variant frequency less than `freq`.
"""
function consensus(refseq::LongDNASeq, vars::AbstractArray{Variant}; freq::Float64=0.5)
    consensus_seq = copy(refseq)

    for i in eachindex(refseq)
        containingvars = filter(v -> varposition(v) == i, vars)
        sort!(containingvars; lt=(x, y) -> alternatedepth(x) < alternatedepth(y))

        if length(containingvars) < 1
            continue
        end #if

        mostvar = last(containingvars)

        if (alternatedepth(mostvar) / totaldepth(mostvar)) > freq
            consensus_seq = mutate(consensus_seq, Haplotype([mostvar]))
        end #if
    end #for

    return consensus_seq
end #function
