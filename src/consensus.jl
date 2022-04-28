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

"""
    consensus_variants(vars::AbstractArray{Variant}; freq::Float64=0.5)

Get a list of [`Variant`](@ref)s from `vars` that are the majority base with a frequency
greater than `freq`
"""
function consensus_variants(vars::AbstractArray{Variant}; freq::Float64=0.5)
    con_vars = Variant[]

    # Loop through every genomic region (CHROM:POS)
    for chrom in unique(chromosome.(vars)), pos in unique(varposition.(vars))
        # Get all the variants that match this particular genomic region
        containingvars = filter(v -> chromosome(v) == chrom && varposition(v) == pos, vars)

        # Sort those variants by depth and pick the most occuring variant
        sort!(containingvars; lt=(x, y) -> alternatedepth(x) < alternatedepth(y))
        mostvar = last(containingvars)

        if (alternatedepth(mostvar) / totaldepth(mostvar)) > freq
            push!(con_vars, mostvar)
        end #if
    end #for

    return con_vars
end #function
