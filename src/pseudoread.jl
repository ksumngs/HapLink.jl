struct Pseudoread
    startpos::UInt64
    endpos::UInt64
    read::Variant
end #struct

function Pseudoread(startpos::Integer, endpos::Integer, read::Variant)
    (startpos < 0 || endpos < 0) && error("Positions cannot be less than zero")
    return Pseudoread(UInt64(startpos), UInt64(endpos), read)
end #function

function Pseudoread(query::Union{SAM.Record,BAM.Record}, reference::NucleotideSeq)
    XAM = query isa SAM.Record ? SAM : BAM

    return Pseudoread(
        XAM.position(query), XAM.rightposition(query), Variant(query, reference)
    )
end #function

_startpos(x::Pseudoread) = x.startpos
_endpos(x::Pseudoread) = x.endpos
_read(x::Pseudoread) = x.read

function _pseudoreads(sam::Union{AbstractString,AbstractPath}, consensus::NucleotideSeq)
    XAM = _issam(sam) ? SAM : BAM

    returned_reads = Pseudoread[]

    reader = open(XAM.Reader, sam)
    for r in reader
        push!(returned_reads, Pseudoread(r, consensus))
    end #for
    close(reader)

    return returned_reads
end #function
