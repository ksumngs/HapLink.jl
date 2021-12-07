export findsimulatedhaplotypes

function findsimulatedhaplotypes(
    variants::AbstractVector{Variant},
    bamfile::AbstractString,
    D::Int,
    α::Float64;
    iterations=1000
)

    # TODO: Abstract the simulated haplotype finding into a higher-level function

    # Find every possible pair of variants. These may be valid haplotypes in their own
    # right, but for right now, we are just going to use them to find linkage between pairs
    variantpairs = combinations(variants, 2)

    # Create a place to store linked pairs and their statistics
    linkedvariantpairhaplotypes = Dict{Haplotype, Matrix{Int}}()

    # Calculate the linkage between every possible variant pair, saving the pair as a
    # haplotype if it exhibited linkage
    for variantpair in variantpairs
        pairedhaplotype = Haplotype(variantpair)
        hapcount = findsimulatedoccurrences(pairedhaplotype, bamfile, iterations=iterations)
        if linkage(hapcount)[2] <= α && last(hapcount) >= D
            linkedvariantpairhaplotypes[pairedhaplotype] = hapcount
        end #if
    end #for

    # Use (abuse?) Julia's awesome broadcasting and spatting to get a non-redundant list of
    # every variant position that exhibited linkage with any other variant position
    linkedvariants = unique(
        cat(map(h -> h.mutations, collect(keys(linkedvariantpairhaplotypes)))..., dims=1)
    )

    # Create a new dict with the basic structure
    # variant => [all possibly linked variants]
    possiblelinkages = Dict{Variant, AbstractArray{Variant}}()

    # Based on the information of what variants exhibit linkage disequilibrium while in
    # pairs, for a variant foo, find every other variant that exhibited paired linkage with
    # foo, and store those other variants in a dict, once again abusing Julia's broadcasting
    # and splatting
    for variant in linkedvariants
        possiblelinkages[variant] = sort(
            unique(
                cat(
                    map(
                        h -> h.mutations,
                        filter(
                            h -> variant in h.mutations,
                            collect(keys(linkedvariantpairhaplotypes))
                        )
                    )...,
                dims=1
                )
            )
        )
    end #for

    # Convert the possible haplotype pairings into non-redundant Haplotype objects
    allvariantcombos = Haplotype.(unique(sort.(values(possiblelinkages))))

    # Set aside a place to put haplotypes that we'll return to the caller
    returnedhaplotypes = Dict{Haplotype, Any}()

    # Calculate the linkage between any new possible haplotypes
    for haplotype in allvariantcombos
        if haskey(linkedvariantpairhaplotypes, haplotype)
            returnedhaplotypes[haplotype] = linkedvariantpairhaplotypes[haplotype]
        else
            hapcount = findsimulatedoccurrences(haplotype, bamfile, iterations=iterations)
            if linkage(hapcount)[2] <= α && last(hapcount) >= D
                returnedhaplotypes[haplotype] = hapcount
            end #if
        end #if
    end #for

    # TODO: Add the single-variant haplotypes back in

    return returnedhaplotypes

end #function
