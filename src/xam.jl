function _xam_record_switch(r::Union{Type{SAM.Record},Type{BAM.Record}})
    return r <: SAM.Record ? :SAM : :BAM
end

function _xam_reader_switch(r::Union{Type{S},Type{B}}) where {S<:SAM.Reader,B<:BAM.Reader}
    return r <: SAM.Reader ? :SAM : :BAM
end #function

function _issam(file::Union{AbstractPath,AbstractString})
    for line in readlines(string(file))
        if startswith(line, '@')
            if startswith(line, "@SQ")
                return true
            end #if
        end #if
    end #for

    return false
end #function

@generated function interval(r::Union{SAM.Record,BAM.Record})
    XAM = _xam_record_switch(r)
    quote
        return Interval($XAM.refname(r), $XAM.position(r), $XAM.rightposition(r))
    end #quote
end #function
