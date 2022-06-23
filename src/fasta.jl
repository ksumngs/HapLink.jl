function _first_record(fasta::Union{AbstractString,AbstractPath})
    reader = FASTA.Reader(open(string(fasta), "r"))
    record = first(reader)
    if !isempty(reader)
        @warn "Reference $fasta contains more than one sequence. Only the first will be used"
    end #if
    close(reader)
    return record
end #function
