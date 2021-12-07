export findsimulatedhaplotypes
export findsimulatedoccurrences

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

"""
    function findsimulatedoccurrences(args...; kwargs...)

Find the number of times a particular haplotype is supported by a maximum likelihood
simulation of combining aligned reads in a BAM file.

# Arguments
- `haplotype::Haplotype`: The combination of variants to test the aligned reads for evidence
    of
- `bamfile::AbstractString`: The path to a BAM file containing aligned reads to be tested
    for instances of `haplotype`

# Keywords
- `iterations::Integer=1000`: The number of times to combine reads and test for the presence
    of `haplotype`

`findsimulatedoccurrences` examines each variant position in `haplotype` and finds a random
read containing that position from `bamfile`. It will then check if the next variant
position is contained on the previous read, and if not pick a new random read that contains
that variant position. In this way, it assembles a set of reads that conceivably could have
come from the same template strand via maximum likelihood.

From that set of reads, `findsimulatedoccurrences` returns an ``N``-dimensional matrix where
``N`` is the number of variant positions in `haplotypes`. The ``1`` index position in the
``n``th dimension represents the number of times the ``n``th variant position was found to
have the reference base called, while the ``2`` index position represents the number of
times the ``n``th variant position was found to have the alternate base called. E.g.
`first(findsimulatedoccurrences(...))` gives the number of times the all-reference base
haplotype was found in the simulation, while `findsimulatedoccurrences(...)[end]` gives the
number of times the all-alternate base haplotype was found.
"""
function findsimulatedoccurrences(
    haplotype::Haplotype,
    bamfile::AbstractString;
    iterations=1000
)

    # Extract the SNPs we care about
    mutations = haplotype.mutations

    # Create an empty array for the simulated long reads
    pseudoreads = Array{Symbol}(undef, iterations, length(mutations))

    # Start reading the BAM file
    open(BAM.Reader, bamfile) do bamreader
        # Collect the reads
        reads = collect(bamreader)

        # Start iterating
        Threads.@threads for i ∈ 1:iterations
            # Get the reads that contain the first mutation
            lastcontainingreads = filter(
                b -> BAM.position(b) < mutations[1].position && BAM.rightposition(b) > mutations[1].position,
                reads
            )

            # Pull a random read from that pool
            lastread = rand(lastcontainingreads)

            # Find this read's basecall at that position
            basecall = baseatreferenceposition(lastread, mutations[1].position)
            basematch = matchvariant(basecall, mutations[1])

            pseudoreads[i, 1] = basematch

            for j ∈ 2:length(mutations)
                if (BAM.position(lastread) < mutations[j].position && BAM.rightposition(lastread) > mutations[j].position)
                    thisread = lastread
                else
                    thiscontainingreads = filter(
                        b -> BAM.position(b) > BAM.rightposition(lastread) && BAM.position(b) < mutations[j].position && BAM.rightposition(b) > mutations[j].position,
                        reads
                    )
                    if length(thiscontainingreads) < 1
                        pseudoreads[i,j] = :other
                        continue
                    end #if
                    thisread = rand(thiscontainingreads)
                end #if

                # Find this read's basecall at that position
                basecall = baseatreferenceposition(thisread, mutations[j].position)
                basematch = matchvariant(basecall, mutations[j])

                pseudoreads[i, j] = basematch

                lastread = thisread
            end #for
        end #for
    end #do

    # Set up haplotype counts
    hapcounts = zeros(Int, repeat([2], length(mutations))...)

    for i ∈ 1:iterations
        matches = pseudoreads[i, :]
        if !any(matches .== :other)
            coordinate = CartesianIndex((Int.(matches .== :alternate) .+ 1)...)
            hapcounts[coordinate] += 1
        end #if
    end #for

    return hapcounts
end #function
