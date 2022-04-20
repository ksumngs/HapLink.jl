"""
    doescontain(snp::SNP, rec::Union{SAM.Record,BAM.Record})

Determines if the location and base mutation of `snp` are present in the sequence of `rec`

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> samrecord = SAM.Record(HapLink.Examples.MutStrings[1]);

julia> doescontain(SNP("ref", 12, DNA_T, DNA_C), samrecord)
false

julia> doescontain(SNP("ref", 17, DNA_T, DNA_A), samrecord)
false

julia> doescontain(SNP("ref", 17, DNA_T, DNA_C), samrecord)
true
```
"""
function doescontain(snp::SNP, rec::Union{SAM.Record,BAM.Record})
    loc = location(snp)

    if !doescontain(loc, rec)
        return false
    end #if

    if first(basesat(loc, rec)) != altbase(snp)
        return false
    end #if

    return true
end #function

"""
    depth(snp::SNP, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Calculate the number of times `snp` is represented within `reads`

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> samrecords = SAM.Record.(HapLink.Examples.MutStrings);

julia> depth(SNP("ref", 17, DNA_T, DNA_C), samrecords)
2
```
"""
function depth(snp::SNP, reads::AbstractVector{T}) where {T<:Union{SAM.Record,BAM.Record}}
    return ThreadsX.count(r -> doescontain(snp, r), reads)
end #function

"""
    mean_quality(snp::SNP, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Calculates the mean PHRED quality of `snp` in the sequences of `reads`.

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> mutrecords = SAM.Record.(HapLink.Examples.MutStrings);

julia> mean_quality(SNP("ref", 12, DNA_T, DNA_T), mutrecords) # reference
30.0

julia> mean_quality(SNP("ref", 17, DNA_T, DNA_C), mutrecords)
35.0
```
"""
function mean_quality(
    snp::SNP, reads::AbstractVector{T}
) where {T<:Union{SAM.Record,BAM.Record}}
    containingreads = filter(r -> doescontain(snp, r), reads)
    loc = location(snp)
    return mean(ThreadsX.map(r -> base_quality(loc, r), containingreads))
end #function

"""
    mean_fractional_position(snp::SNP, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Finds the mean position as a fraction  between 0 and 1 of `snp` within a set of `reads`.

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> mean_fractional_position(SNP("ref", 9, DNA_A, DNA_A), SAM.Record.(HapLink.Examples.MutStrings))
0.3358798064680418

julia> mean_fractional_position(SNP("ref", 17, DNA_T, DNA_C), SAM.Record.(HapLink.Examples.MutStrings))
0.4732620320855615
```
"""
function mean_fractional_position(
    snp::SNP, reads::AbstractVector{T}
) where {T<:Union{SAM.Record,BAM.Record}}
    containingreads = filter(r -> doescontain(snp, r), reads)
    loc = location(snp)
    return mean(ThreadsX.map(r -> fractional_position(loc, r), containingreads))
end #function

"""
    frequency(snp::SNP, reads::AbstractVector{T}) where T <: Union{SAM.Record,BAM.Record}

Calculate the relative depth of `snp` within `reads`.

# Example

```jldoctest
julia> using BioSymbols, GenomicFeatures, XAM

julia> frequency(SNP("ref", 17, DNA_T, DNA_C), SAM.Record.(HapLink.Examples.MutStrings))
0.6666666666666666
```
"""
function frequency(
    snp::SNP, reads::AbstractVector{T}
) where {T<:Union{SAM.Record,BAM.Record}}
    return depth(snp, reads) / depth(location(snp), reads)
end #function
