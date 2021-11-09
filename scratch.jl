using Combinatorics

function islinked(combo::String)
    # Arbitrary linkage simulator
    # A-B-D
    # C-E
    # F
    # G-H-I-J
    linkages = ["A-B-D", "C-E", "F-K", "G-H-I-J", "K-L", "M"]
    alleles = string.(split(combo, "-"))

    for link in linkages
        if all(occursin.(alleles, link))
            return true
        end
    end
    return false
end #function

variants = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]

varcombos = combinations(variants, 2)

linkedvarcombos = []

for varcombo in varcombos
    var = join(varcombo, "-")
    if islinked(var)
        push!(linkedvarcombos, varcombo)
    end
end

linkedvars = unique(cat(linkedvarcombos..., dims = 1))

varhaps = Dict()

for l in linkedvars
    varhaps[l] = sort(unique(cat(filter(h -> l in h, linkedvarcombos)..., dims=1)))
end

newvarcombos = unique(values(varhaps))

filter!(v -> islinked(join(v, "-")), newvarcombos)
