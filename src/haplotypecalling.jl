export find_haplotypes
export longread_genome
export simulate_genome
export occurrence_matrix
export linkage
export sumsliced

"""
    find_haplotypes(
        variants::AbstractVector{Variant},
        bamfile::AbstractString,
        D::Int,
        α::Float64,
        haplotypemethod
    )

Find all combinations of `variants` that the reads in `bamfile` support as a valid haplotype
with a minimum depth of `D` and Χ-squared linkage disequilibrium significance of `α`, with
the haplotypes being converted into genomes via `haplotypemethod`.

# Arguments
- `variants::AbstractVector{Variant}`: A `Vector` of `Variant` objects which will be
    combined into `Haplotype` objects and tested for validity
- `bamfile::AbstractString`: The path to a BAM file containing reads to check for the
    presence of haplotypes
- `D::Int`: The minimum number of times a haplotype must be present in the reads according
    to `haplotypemethod`
- `α::Float64`: The maximum ``Χ``-squared ``p``-value at which to consider a haplotype
    significant and worth returning.
- `haplotypemethod`: A function handle with the signature `f(h::Haplotype,
    b::AbstractString)` which with return a table of the variant basecall matches found.
    Both [`longread_genome`](@ref) and [`simulate_genome`](@ref) fulfil this requirement.

# Returns
- `Dict{Haplotype,Matrix{Int}}`: A dictionary containing every significant haplotype and its
    incidence matrix
"""
function find_haplotypes(
    variants::AbstractVector{Variant},
    bamfile::AbstractString,
    D::Int,
    α::Float64,
    haplotypemethod,
)

    # Find every possible pair of variants. These may be valid haplotypes in their own
    # right, but for right now, we are just going to use them to find linkage between pairs
    variantpairs = combinations(variants, 2)

    # Create a place to store linked pairs and their statistics
    linkedvariantpairhaplotypes = Dict{Haplotype,Matrix{Int}}()

    # Calculate the linkage between every possible variant pair, saving the pair as a
    # haplotype if it exhibited linkage
    for variantpair in variantpairs
        pairedhaplotype = Haplotype(variantpair)
        hapcount = occurrence_matrix(haplotypemethod(pairedhaplotype, bamfile))
        if linkage(hapcount)[2] <= α && last(hapcount) >= D
            linkedvariantpairhaplotypes[pairedhaplotype] = hapcount
        end #if
    end #for

    # Use (abuse?) Julia's awesome broadcasting and spatting to get a non-redundant list of
    # every variant position that exhibited linkage with any other variant position
    linkedvariants = unique(
        cat(map(h -> h.mutations, collect(keys(linkedvariantpairhaplotypes)))...; dims=1)
    )

    # Create a new dict with the basic structure
    # variant => [all possibly linked variants]
    possiblelinkages = Dict{Variant,AbstractArray{Variant}}()

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
                            collect(keys(linkedvariantpairhaplotypes)),
                        ),
                    )...;
                    dims=1,
                ),
            ),
        )
    end #for

    # Convert the possible haplotype pairings into non-redundant Haplotype objects
    allvariantcombos = Haplotype.(unique(sort.(values(possiblelinkages))))

    # Set aside a place to put haplotypes that we'll return to the caller
    returnedhaplotypes = Dict{Haplotype,Any}()

    # Calculate the linkage between any new possible haplotypes
    for haplotype in allvariantcombos
        if haskey(linkedvariantpairhaplotypes, haplotype)
            returnedhaplotypes[haplotype] = linkedvariantpairhaplotypes[haplotype]
        else
            hapcount = occurrence_matrix(haplotypemethod(haplotype, bamfile))
            if linkage(hapcount)[2] <= α && last(hapcount) >= D
                returnedhaplotypes[haplotype] = hapcount
            end #if
        end #if
    end #for

    # Add the single-variant haplotypes back in
    confirmedlinkedvariants = unique(
        cat(mutations.(collect(keys(returnedhaplotypes)))...; dims=1)
    )
    singlevariants = filter(v -> !(v in confirmedlinkedvariants), variants)
    for var in singlevariants
        varhap = Haplotype(var)
        altdepth = alternatedepth(var)
        refdepth = totaldepth(var) - altdepth
        returnedhaplotypes[varhap] = [refdepth altdepth]
    end #for

    return returnedhaplotypes
end #function

"""
    longread_genome(haplotype::Haplotype, bamfile::AbstractString)

Parse the whole-genome length reads in `bamfile` to determine each read's basecall at every
variant position within `haplotype`.

# Arguments
- `haplotype::Haplotype`: The combination of variants to test basecalls against
- `bamfile::AbstractString`: The path to a BAM file containing aligned reads to be tested
    for instances of `haplotype`

# Returns
- `MxN Array{Symbol}` where `M` is the number of reads present in `bamfile` and
    `N=length(haplotype.mutations)`: A table of which base each read has in every variant
    position of `haplotype`. The table has reads for rows and variant positions for columns,
    e.g. longread_genome(...)[5,2] gives the basecall for the fifth read at the second
    variant position. Basecalls are given as `Symbol` objects with possible values of
    - `:reference`
    - `:alternate`
    - `:other`
"""
function longread_genome(haplotype::Haplotype, bamfile::AbstractString)

    # Extract the SNPs we care about
    mutations = haplotype.mutations

    # Create the return object
    pseudoreads = Symbol[]

    # Start reading the BAM file
    open(BAM.Reader, bamfile) do bamreader
        # Collect the reads
        reads = collect(bamreader)

        # Get only the reads that contain all of the variant positions
        containingreads = filter(
            b ->
                BAM.position(b) < min(varposition.(mutations)...) &&
                    BAM.rightposition(b) > max(varposition.(mutations)...),
            reads,
        )

        # Piece out if there aren't any full-length reads
        if length(containingreads) < 1
            @warn "No reads in $bamfile contained all variant positions. Exiting now."
            return reshape(repeat([:other], length(mutations)), (1, length(mutations)))
        end #if

        # Overwrite the return object as the correct dimensions
        pseudoreads = Array{Symbol}(undef, length(containingreads), length(mutations))

        # Check every NGS read that contains all positions
        for i in 1:length(containingreads)
            # Make it like a foreach loop
            record = containingreads[i]

            # Check agains every variant position
            for j in 1:length(mutations)
                # Make it like a foreach loop
                mutation = mutations[j]

                # Pull the basecall
                basecall = baseatreferenceposition(record, varposition(mutation))
                basematch = matchvariant(basecall, mutation)

                # Put into the results
                pseudoreads[i, j] = basematch
            end #for
        end #for
    end #do

    return pseudoreads
end #function

"""
    simulate_genome(haplotype::Haplotype, bamfile::AbstractString; iterations::Int64=1000)

Simulate a set of `iterations` reads that contain all of the variant positions in
`haplotype` containing the actual reads present in `bamfile` via a maximum likelihood
method.

`simulate_genome` examines each variant position in `haplotype` and finds a random read
containing that position from `bamfile`. It will then check if the next variant position is
contained on the previous read, and if not pick a new random read that contains that variant
position. In this way, it assembles a set of reads that conceivably could have come from the
same template strand via maximum likelihood.

# Arguments
- `haplotype::Haplotype`: The combination of variants to test the aligned reads for evidence
    of
- `bamfile::AbstractString`: The path to a BAM file containing aligned reads to be tested
    for instances of `haplotype`

# Keywords
- `iterations::Integer=1000`: The number of times to combine reads and test for the presence
    of `haplotype`

# Returns
- `MxN Array{Symbol}` where `M=iterations` and `N=length(haplotype.mutations)`: A table of
    which base each simulated read has in every variant position of `haplotype`. The table
    has reads for rows and variant positions for columns, e.g. simulate_genome(...)[5,2]
    gives the basecall for the fifth simulated read at the second variant position.
    Basecalls are given as `Symbol` objects with possible values of
    - `:reference`
    - `:alternate`
    - `:other`
"""
function simulate_genome(haplotype::Haplotype, bamfile::AbstractString; iterations=1000)

    # TODO: implement an overlapped-read ML algorithm

    # Extract the SNPs we care about
    mutations = haplotype.mutations

    # Create an empty array for the simulated long reads
    pseudoreads = Array{Symbol}(undef, iterations, length(mutations))

    # Start reading the BAM file
    open(BAM.Reader, bamfile) do bamreader
        # Collect the reads
        reads = collect(bamreader)

        # Start iterating
        Threads.@threads for i in 1:iterations
            # Get the reads that contain the first mutation
            lastcontainingreads = filter(
                b ->
                    BAM.position(b) < mutations[1].position &&
                        BAM.rightposition(b) > mutations[1].position,
                reads,
            )

            # Pull a random read from that pool
            lastread = rand(lastcontainingreads)

            # Find this read's basecall at that position
            basecall = baseatreferenceposition(lastread, mutations[1].position)
            basematch = matchvariant(basecall, mutations[1])

            pseudoreads[i, 1] = basematch

            for j in 2:length(mutations)
                if (
                    BAM.position(lastread) < mutations[j].position &&
                    BAM.rightposition(lastread) > mutations[j].position
                )
                    thisread = lastread
                else
                    thiscontainingreads = filter(
                        b ->
                            BAM.position(b) > BAM.rightposition(lastread) &&
                                BAM.position(b) < mutations[j].position &&
                                BAM.rightposition(b) > mutations[j].position,
                        reads,
                    )
                    if length(thiscontainingreads) < 1
                        pseudoreads[i, j] = :other
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

    return pseudoreads
end #function

"""
    occurrence_matrix(readmatches::AbstractArray{Symbol})

Transforms the haplotype occurrence table `readmatches` into an incidence matrix

# Arguments
- `readmatches::AbstractArray{Symbol}`: An ``m``x``n`` array where ``m`` is the number of
    reads represented, and ``n`` is the number of variants in the haplotype considered, e.g.
    readmatches[4,3] represents the match value for the third variant in the fourth read.
    Valid values in the array are `:reference`, `:alternate`, and `:other`.

# Returns
- `2x2x... Array{Int64, size(readmatches)[2]}`: An ``N``-dimensional matrix where ``N`` is
    the number of variant positions in `readmatches`. The ``1`` index position in the
    ``n``th dimension represents the number of times the ``n``th variant position was found
    to have the reference base called, while the ``2`` index position represents the number
    of times the ``n``th variant position was found to have the alternate base called. E.g.
    `first(occurrence_matrix(reads))` gives the number of times the all-reference base
    haplotype was found in `reads`, while `occurrence_matrix(reads)[end]` gives the number
    of times the all-alternate base haplotype was found.

# Example
```jldoctest
julia> pseudoreads = [
           :reference :reference :reference
           :reference :reference :alternate
           :reference :reference :alternate
           :reference :reference :other
       ];

julia> occurrence_matrix(pseudoreads)
2×2×2 Array{Int64, 3}:
[:, :, 1] =
 1  0
 0  0

[:, :, 2] =
 2  0
 0  0
```
"""
function occurrence_matrix(readmatches::AbstractArray{Symbol})
    # Set up an n-dimensional matrix to store haplotype counts
    hapcounts = zeros(Int, repeat([2], size(readmatches)[2])...)

    # Transform each match row into a count in the incidence matrix
    for i in 1:size(readmatches)[1]
        matches = readmatches[i, :]
        if !any(matches .== :other)
            coordinate = CartesianIndex((Int.(matches .== :alternate) .+ 1)...)
            hapcounts[coordinate] += 1
        end #if
    end #for

    return hapcounts
end #function

"""
    linkage(counts::AbstractArray{Int})

Calculates the linkage disequilibrium and Chi-squared significance level of a combination of
haplotypes whose number of occurrences are given by `counts`.

`counts` is an ``N``-dimensional array where the ``N``th dimension represents the ``N``th
variant call position within a haplotype. `findoccurrences` produces such an array.
"""
function linkage(counts::AbstractArray{Int})
    # Get the probability of finding a perfect reference sequence
    P_allref = first(counts) / sum(counts)

    # Get the probabilities of finding reference bases in any of the haplotypes
    P_refs = sumsliced.([counts], 1:ndims(counts)) ./ sum(counts)

    # Calculate linkage disequilibrium
    Δ = P_allref - prod(P_refs)

    # Calculate the test statistic
    r = Δ / (prod(P_refs .* (1 .- P_refs))^(1 / ndims(counts)))
    Χ_squared = r^2 * sum(counts)

    # Calculate the significance
    p = 1 - cdf(Chisq(1), Χ_squared)

    return Δ, p
end #function

"""
    sumsliced(A::AbstractArray, dim::Int, pos::Int=1)

Sum all elements that are that can be referenced by `pos` in the `dim` dimension of `A`.

# Example

```jldoctest
julia> A = reshape(1:8, 2, 2, 2)
2×2×2 reshape(::UnitRange{Int64}, 2, 2, 2) with eltype Int64:
[:, :, 1] =
 1  3
 2  4

[:, :, 2] =
 5  7
 6  8

julia> sumsliced(A, 2)
14

julia> sumsliced(A, 2, 2)
22
```

Heavily inspired by Holy, Tim "Multidimensional algorithms and iteration"
<https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes>
"""
function sumsliced(A::AbstractArray, dim::Int, pos::Int=1)
    i_pre  = CartesianIndices(size(A)[1:(dim - 1)])
    i_post = CartesianIndices(size(A)[(dim + 1):end])
    return sum(A[i_pre, pos, i_post])
end #function

function serialize_yaml(h::Pair; reason::Union{String,Nothing}=nothing)
    occurrences = "occurrences:\n"
    if sum(size(h.second)) > 3
        for i in CartesianIndices(h.second)
            location = [Tuple(i)...]
            variantpattern = string.(replace(replace(location, 1 => "ref"), 2 => "alt"))
            key = join(variantpattern, "_")
            occurrences = string(occurrences, "  ", key, ": ", h.second[i], "\n")
        end #for
    else
        occurrences = string(
            occurrences, "  ref: ", h.second[1], "\n  alt: ", h.second[2], "\n"
        )
    end #if

    return string(
        serialize_yaml(h.first; reason=reason),
        occurrences,
        "Δ: ",
        linkage(h.second)[1],
        "\n",
        "p: ",
        linkage(h.second)[2],
        "\n",
    )
end #function
