@kwdef mutable struct HaplotypeMeta
    frequency::Union{Float64,Nothing} = nothing
    linkage::Union{Float64,Nothing} = nothing
    significance::Union{Float64,Nothing} = nothing
    name::Union{String,Nothing} = nothing
end #struct
