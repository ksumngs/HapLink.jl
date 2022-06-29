"""
    _indexed_reader(bam::Union{AbstractPath,AbstractString}, int::Interval)

Provides an iterator over records in `bam` that overlap with `int`. Note that if no index
file is available for `bam`, then the return value will fall back to an iterator of all
records in `bam`.
"""
function _indexed_reader(bam::Union{AbstractPath,AbstractString}, int::Interval)
    # Find the index file
    bai = _find_bam_index(bam)

    if isnothing(bai)
        # There is no index file
        reader = _is_sam ? SAM.Reader : BAM.Reader
        return (reader, reader(open(string(bam), "r")))
    else
        # There is an index, so we can zip through this
        reader = open(BAM.Reader, string(bam); index=string(bai))
        return (BAM.OverlapIterator, eachoverlap(reader, int))
    end
end
