using GenomicFeatures
using Statistics
using XAM

export basesat
export base_quality
export depth
export doescontain
export fractional_position
export mean_quality

"""
    basesat(int::GenomicFeatures.Interval, rec::Union{SAM.Record, BAM.Record})

Gets the sequence of `rec` at the reference position of `int`.

# Example

```jldoctest
julia> using BioSequences, BioSymbols, GenomicFeatures, XAM

julia> samrecord = SAM.Record(HapLink.Examples.SAMStrings[1]);

julia> location = Interval("ref", 15, 17);

julia> basesat(location, samrecord)
3nt DNA Sequence:
GAT
```
"""
@generated function basesat(int::Interval, rec::Union{BAM.Record,SAM.Record})
    XAM = rec <: SAM.Record ? :SAM : :BAM

    quote
        lpos = leftposition(int)
        rpos = rightposition(int)
        seqsize = rpos - lpos + 1

        alignment = $XAM.alignment(rec)
        sequence = $XAM.sequence(rec)

        outbases = Vector{NucleicAcid}(undef, seqsize)

        for (i, pos) in enumerate(lpos:rpos)
            refpos = first(ref2seq(alignment, pos))
            outbases[i] = sequence[refpos]
        end #for

        NucType = typeof(sequence)
        return NucType(outbases)
    end #quote
end #function

"""
    doescontain(int::GenomicFeatures.Interval, rec::Union{SAM.Record,BAM.Record})

Determines if `int` is fully included within the sequence of `rec` as matching bases.

Named `doescontain` to avoid name conflict with `Base.contains()`

# Example

```jldoctest
julia> using GenomicFeatures, XAM

julia> doescontain(Interval("ref", 10, 17), SAM.Record(HapLink.Examples.SAMStrings[1]))
true

julia> doescontain(Interval("ref", 14, 25), SAM.Record(HapLink.Examples.SAMStrings[1]))
false

julia> doescontain(Interval("ref", 30, 31), SAM.Record(HapLink.Examples.SAMStrings[1]))
false

julia> doescontain(Interval("ref", 23, 26), SAM.Record(HapLink.Examples.SAMStrings[4]))
false
```
"""
@generated function doescontain(int::Interval, rec::Union{SAM.Record,BAM.Record})
    XAM = rec <: SAM.Record ? :SAM : :BAM

    quote
        if seqname(int) != $XAM.refname(rec)
            return false
        end #if

        lpos = leftposition(int)
        if lpos < $XAM.position(rec)
            return false
        end #if

        rpos = rightposition(int)
        if rightposition(int) > $XAM.rightposition(rec)
            return false
        end #if

        aln = $XAM.alignment(rec)
        for pos in lpos:rpos
            if !ismatchop(last(ref2seq(aln, pos)))
                return false
            end #if
        end #for

        return true
    end #quote
end #function

"""
    depth(int::Interval, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Calculate the number of `Record`s in `reads` that contain `int`

# Example

```jldoctest
julia> using GenomicFeatures, XAM

julia> samrecords = SAM.Record.(HapLink.Examples.SAMStrings);

julia> depth(Interval("ref", 17, 18), samrecords)
3
```
"""
function depth(
    int::Interval, reads::AbstractVector{T}
) where {T<:Union{SAM.Record,BAM.Record}}
    return count(r -> doescontain(int, r), reads)
end #function

"""
    base_quality(int::Interval, rec::Union{SAM.Record,BAM.Record})

Calculates the mean PHRED quality across all of the basecalls in `int` for the read in
`rec`. Only considers matching alignments, i.e. insertions and deletions are removed

# Example

```jldoctest
julia> using GenomicFeatures, XAM

julia> samrecord = SAM.Record(HapLink.Examples.SAMStrings[1]);

julia> base_quality(Interval("ref", 17, 17), samrecord)
30.0

julia> base_quality(Interval("ref", 7, 21), samrecord)
30.0
```
"""
@generated function base_quality(int::Interval, rec::Union{SAM.Record,BAM.Record})
    XAM = rec <: SAM.Record ? :SAM : :BAM

    quote
        lpos = leftposition(int)
        rpos = rightposition(int)
        seqsize = rpos - lpos + 1

        alignment = $XAM.alignment(rec)
        qualities = $XAM.quality(rec)

        outquals = Vector{Union{Missing,Float64}}(missing, seqsize)

        for (i, pos) in enumerate(lpos:rpos)
            refpos, op = ref2seq(alignment, pos)
            if ismatchop(op)
                outquals[i] = qualities[refpos]
            end #if
        end #for

        return mean(collect(skipmissing(outquals)))
    end #quote
end #function

"""
    mean_quality(int::Interval, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Calculates the mean PHRED quality for `int` in the sequences of `reads`. Passing `int`s of
more than one position will average the per-base scores and the per-read scores, and may
yield unexpected results. Records within `reads` that do not contain `int` will be ignored.

# Example

```jldoctest
julia> using GenomicFeatures, XAM

julia> mutrecords = SAM.Record.(HapLink.Examples.MutStrings);

julia> mean_quality(Interval("ref", 12, 12), mutrecords)
30.0

julia> mean_quality(Interval("ref", 17, 17), mutrecords)
33.333333333333336
```
"""
function mean_quality(
    int::Interval, reads::AbstractVector{T}
) where {T<:Union{SAM.Record,BAM.Record}}
    containingreads = filter(r -> doescontain(int, r), reads)
    return mean(base_quality.([int], containingreads))
end #function

"""
    fractional_position(int::Interval, rec::Union{SAM.Record,BAM.Record})

Get the position of `int` as a fraction within the sequence of `rec`.

# Example

```jldoctest
julia> using GenomicFeatures, XAM

julia> samrecord = SAM.Record.(HapLink.Examples.SAMStrings[5]);

julia> fractional_position(Interval("ref", 29, 29), samrecord)
0.2

julia> fractional_position(Interval("ref", 29, 33), samrecord)
0.6

julia> fractional_position(Interval("ref", 33, 33), samrecord)
1.0
```
"""
function fractional_position(int::Interval, rec::SAM.Record)
    leftpos = first(ref2seq(SAM.alignment(rec), leftposition(int)))
    rightpos = first(ref2seq(SAM.alignment(rec), rightposition(int)))
    avgpos = mean([leftpos, rightpos])
    return avgpos / SAM.seqlength(rec)
end #function

function fractional_position(int::Interval, rec::BAM.Record)
    leftpos = first(ref2seq(BAM.alignment(rec), leftposition(int)))
    rightpos = first(ref2seq(BAM.alignment(rec), rightposition(int)))
    avgpos = mean([leftpos, rightpos])
    return avgpos / BAM.seqlength(rec)
end #function

"""
    fractional_position(int::Interval, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Finds the mean position as a fraction between 0 and 1 of `int` within a set of `reads`.
Reads without `int` are ignored.

# Example

```jldoctest
julia> using GenomicFeatures, XAM

julia> fractional_position(Interval("ref", 9, 9), SAM.Record(HapLink.Examples.SAMStrings[1]))
0.17647058823529413

julia> fractional_position(Interval("ref", 9, 9), SAM.Record(HapLink.Examples.SAMStrings[2]))
0.2857142857142857

julia> fractional_position(Interval("ref", 9, 9), SAM.Record(HapLink.Examples.SAMStrings[3]))
0.5454545454545454

julia> fractional_position(Interval("ref", 9, 9), SAM.Record.(HapLink.Examples.SAMStrings))
0.3358798064680418
```
"""
function fractional_position(
    int::Interval, reads::AbstractVector{T}
) where {T<:Union{SAM.Record,BAM.Record}}
    containingreads = filter(r -> doescontain(int, r), reads)
    return mean(fractional_position.([int], containingreads))
end #function
