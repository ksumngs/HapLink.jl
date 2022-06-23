function _xam_record_switch(r::Union{Type{SAM.Record},Type{BAM.Record}})
    return r <: SAM.Record ? :SAM : :BAM
end
