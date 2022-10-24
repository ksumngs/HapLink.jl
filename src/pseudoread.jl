struct Pseudoread
    startpos::UInt64
    endpos::UInt64
    read::Variant
end #struct

function Pseudoread(startpos::Integer, endpos::Integer, read::Variant)
    (startpos < 0 || endpos < 0) && error("Positions cannot be less than zero")
    return Pseudoread(UInt64(startpos), UInt64(endpos), read)
end #function

_startpos(x::Pseudoread) = x.startpos
_endpos(x::Pseudoread) = x.endpos
_read(x::Pseudoread) = x.read
