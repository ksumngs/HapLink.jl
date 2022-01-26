using GenomicFeatures
using XAM

export basesat
export doescontain

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

julia> # Create a record: this example comes from the SAMv1 spec

julia> samrecord = SAM.Record("r001\\t99\\tref\\t7\\t30\\t8M2I4M1D3M\\t=\\t37\\t39\\tTTAGATAAAGGATACTG\\t*")
XAM.SAM.Record:
    template name: r001
             flag: 99
        reference: ref
         position: 7
  mapping quality: 30
            CIGAR: 8M2I4M1D3M
   next reference: =
    next position: 37
  template length: 39
         sequence: TTAGATAAAGGATACTG
     base quality: <missing>
   auxiliary data:

julia> location = Interval("ref", 15, 17);

julia> basesat(location, samrecord)
3nt DNA Sequence:
GAT
```
"""
function basesat(int::Interval, rec::SAM.Record)
    poss = myref2seq.([SAM.alignment(rec)], leftposition(int):rightposition(int))
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

julia> # Use the SAM spec example record

julia> samrecord = SAM.Record("r001\\t99\\tref\\t7\\t30\\t8M2I4M1D3M\\t=\\t37\\t39\\tTTAGATAAAGGATACTG\\t*");

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
