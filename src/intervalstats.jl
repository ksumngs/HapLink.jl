using GenomicFeatures
using Statistics
using XAM

export basesat
export depth
export doescontain
export fractional_position
export mean_quality

#=
samstrings = [
    "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*",
    "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t*",
    "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t*\tSA:Z:ref,29,-,6H5M,17,0;",
    "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCTTCAGC\t*",
    "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t*\tSA:Z:ref,9,+,5S6M,30,1;",
    "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t*\tNM:i:1"
]

refstring = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT"
=#

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
function basesat(int::Interval, rec::SAM.Record)
    poss = ref2seq.([SAM.alignment(rec)], leftposition(int):rightposition(int))
    nuctype = typeof(SAM.sequence(rec))
    return nuctype(map(i -> SAM.sequence(rec)[i], first.(poss)))
end #function

function basesat(int::Interval, rec::BAM.Record)
    poss = myref2seq.([BAM.alignment(rec)], leftposition(int):rightposition(int))
    nuctype = typeof(BAM.sequence(rec))
    return nuctype(map(i -> BAM.sequence(rec)[i], first.(poss)))
end #function

"""
    doescontain(int::GenomicFeatures.Interval, rec::Union{SAM.Record,BAM.Record})

Determines if `int` is fully included within the sequence of `rec`.

Named `doescontain` to avoid name conflict with `Base.contains()`

# Example

This alignment is taken from the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

```plaintext
Coor    12345678901234  5678901234567890123456789012345
ref     AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT

+r001/1       TTAGATAAAGGATACTG
```

| Interval  | Overlap | Result  |
| --------- | ------- | ------- |
| ref:10-17 | Full    | `true`  |
| ref:14-25 | Partial | `false` |
| ref:30-31 | None    | `false` |

```jldoctest
julia> using GenomicFeatures, XAM

julia> samrecord = SAM.Record(HapLink.Examples.SAMStrings[1]);

julia> doescontain(Interval("ref", 10, 17), samrecord)
true

julia> doescontain(Interval("ref", 14, 25), samrecord)
false

julia> doescontain(Interval("ref", 30, 31), samrecord)
false
```
"""
function doescontain(int::Interval, rec::SAM.Record)
    return seqname(int) == SAM.refname(rec) &&
           leftposition(int) >= SAM.position(rec) &&
           rightposition(int) <= SAM.rightposition(rec)
end #function

function doescontain(int::Interval, rec::BAM.Record)
    return seqname(int) == BAM.refname(rec) &&
           leftposition(int) >= BAM.position(rec) &&
           rightposition(int) <= BAM.rightposition(rec)
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
function depth(int::Interval, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}
    return count(r -> doescontain(int, r), reads)
end #function

"""
    mean_quality(int::Interval, rec::Union{SAM.Record,BAM.Record})

Calculates the mean PHRED quality across all of the basecalls in `int` for the read in
`rec`. Only considers matching alignments, i.e. insertions and deletions are removed

# Example

```jldoctest
julia> using GenomicFeatures, XAM

julia> samrecord = SAM.Record(HapLink.Examples.SAMStrings[1]);

julia> mean_quality(Interval("ref", 17, 17), samrecord)
30.0

julia> mean_quality(Interval("ref", 7, 21), samrecord)
30.0
```
"""
function mean_quality(int::Interval, rec::SAM.Record)
    poss = ref2seq.([SAM.alignment(rec)], leftposition(int):rightposition(int))
    matchposs = filter(p -> ismatchop(last(p)), poss)
    return mean(SAM.quality(rec)[first.(matchposs)])
end #function

function mean_quality(int::Interval, rec::BAM.Record)
    poss = ref2seq.([BAM.alignment(rec)], leftposition(int):rightposition(int))
    matchposs = filter(p -> ismatchop(last(p)), poss)
    return mean(BAM.quality(rec)[first.(matchposs)])
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
function mean_quality(int::Interval, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}
    containingreads = filter(r -> doescontain(int, r), reads)
    return mean(mean_quality.([int], containingreads))
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
