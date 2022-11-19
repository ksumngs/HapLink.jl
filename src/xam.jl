function _XAM_(r::Union{SAM.Record,BAM.Record,SAM.Reader,BAM.Reader})
    if r isa SAM.Record || r isa SAM.Reader
        return SAM
    elseif r isa BAM.Record || r isa BAM.Reader
        return BAM
    else
        error("Type checking failed for $(r)")
    end #if
end #function

function _XAM_(f::Union{AbstractString,AbstractPath})
    return _issam(f) ? SAM : BAM
end #function

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

FilePaths.@compat function _find_bam_index(bam::AbstractPath)
    # Check for a '.bam.bai' extension
    bai = Path("$bam.bai")
    isfile(bai) && return bai

    # Check for a swapped extension
    bai = Path("$(first(splittext(bam))).bai")
    isfile(bai) && return bai

    # Couldn't fine an index
    @warn "Couldn't find an index file for $bam. Analysis may be significantly slower."
    return nothing
end #function
