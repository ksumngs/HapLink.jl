Base.@kwdef mutable struct HaplotypeMeta
    frequency::Union{Float64,Nothing} = nothing
    linkage::Union{Float64,Nothing} = nothing
    significance::Union{Float64,Nothing} = nothing
    name::Union{String,Nothing} = nothing
end #struct

function Dict(hp::Pair{Haplotype,HaplotypeMeta})
    return Dict(
        "frequency" => hp[2].frequency,
        "linkage" => hp[2].linkage,
        "significance" => hp[2].significance,
        "name" => hp[2].name,
        "snps" => Dict.(hp[1].mutations),
    )
end #function
