function _xam_record_switch(r::Union{Type{SAM.Record},Type{BAM.Record}})
    return r <: SAM.Record ? :SAM : :BAM
end

function _xam_reader_switch(r::Union{Type{SAM.Reader},Type{BAM.Reader}})
    return r <: SAM.Reader ? :SAM : :BAM
end #function
