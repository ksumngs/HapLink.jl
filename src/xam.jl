_xam_switch(r::Union{Type{SAM.Record},Type{BAM.Record}}) = r <: SAM.Record ? :SAM : :BAM
