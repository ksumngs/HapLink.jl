"""
    findset(lst::AbstractArray{T}, comparator::Function) where {T}

Finds every possible set of items in `lst` where `comparator` returns true for the set.

# Example
```jldoctest
function divisible_by_same_number(x, i)
    for j in 2:i
        if all(y -> y % j == 0, x)
            return true
        end
    end
    return false
end
findset(1:12, x -> divisible_by_same_number(x, 5))

# output
6-element Vector{Vector{Int64}}:
 [2, 4, 6, 8, 10, 12]
 [6, 3, 9, 12]
 [5, 10]
 [1]
 [7]
 [11]
```
"""
function findset(lst::AbstractArray{T}, comparator::Function) where {T}
    # Get all posible pairings of items
    pairs = combinations(lst, 2)

    # Create a place to store valid pairings
    valid_pairs = Vector{Vector{T}}()

    # Search for valid pairings
    for pair in pairs
        if comparator(pair)
            push!(valid_pairs, pair)
        end #if
    end #for

    # Check if there are no valid pairs and abort
    isempty(valid_pairs) && return [[i] for i in lst]

    # Get every item that has a pair match
    set_items = unique(reduce(vcat, valid_pairs))

    # Create a place to store valid sets larger than pairs
    valid_sets = Vector{Vector{T}}()

    # Find all possible sets from every item that was paired
    for paired_item in set_items
        # Find potential matches based on previous pairing results
        matched_pairs = filter(
            x -> first(x) == paired_item || last(x) == paired_item, valid_pairs
        )
        matched_items = unique(reduce(vcat, matched_pairs))

        # Abort if all sets are simply pairs
        length(matched_items) > 2 || continue

        # Starting from the largest possible sets, check every possible combination
        for n in length(matched_items):-1:3
            possible_sets = combinations(matched_items, n)

            for possible_set in possible_sets
                # This set has already been verified as part of a larger set, so don't
                # recalculate
                any(s -> issubset(possible_set, s), valid_sets) && continue

                # Compare and test if this is a valid set
                if comparator(possible_set)
                    push!(valid_sets, possible_set)
                end #if
            end #for
        end #for
    end #for

    # Add any pairs that have not been supersceded to the return set
    for pair in valid_pairs
        any(s -> issubset(pair, s), valid_sets) && continue
        push!(valid_sets, pair)
    end #for

    # Add all lone items in
    lone_items = filter(s -> !(s in set_items), lst)
    for t in lone_items
        push!(valid_sets, [t])
    end #for

    return valid_sets
end #function
