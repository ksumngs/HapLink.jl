using DataFrames

export countbasestats
export transformbamcounts

"""
    countbasestats(bamfile::String, reffile::String)

Count and calculate statistics on the basecalls of the alignment in `bamfile` to the
reference genome in `reffile`. Returns a `DataFrame` with stats on every base in every
alignment position. See [`transformbamcounts`](@ref) for a complete description of the
output DataFrame schema.
"""
function countbasestats(bamfile::String, reffile::String)
    bamanalysis = ""
    open(FASTA.Reader, reffile) do refreader
        for refrecord in refreader
            chromosome = FASTA.identifier(refrecord)
            seqlength = FASTA.seqlen(refrecord)
            bamanalysis = string(
                bamanalysis,
                readchomp(`bam-readcount -f $reffile $bamfile "$chromosome:1-$seqlength"`),
            )
        end #for
    end #do
    return transformbamcounts(string.(split(bamanalysis, "\n")))
end #function

"""
    transformbamcounts(bamcounts::AbstractVector{String})

Convert the output from [bam-readcount](https://github.com/genome/bam-readcount) to a
`DataFrame`.

### Schema

| Column name                            | Description                                                                                    |
| -------------------------------------: | :--------------------------------------------------------------------------------------------- |
| `chr`                                  | Chromosome                                                                                     |
| `position`                             | Position                                                                                       |
| `reference_base`                       | Reference base                                                                                 |
| `depth`                                | Total read depth at `position`                                                                 |
| `base`                                 | Alternate base                                                                                 |
| `count`                                | Number of reads containing `base` at `position`                                                |
| `avg_mapping_quality`                  | Mean mapping quality                                                                           |
| `avg_basequality`                      | Mean base quality at `position`                                                                |
| `avg_se_mapping_quality`               | Mean single-ended mapping quality                                                              |
| `num_plus_strand`                      | Number of reads on the forward strand (N/A)                                                    |
| `num_minus_strand`                     | Number of reads on the reverse strand (N/A)                                                    |
| `avg_pos_as_fraction`                  | Average position on the read as a fraction, 0 (end) to 1 (center)                              |
| `avg_num_mismatches_as_fraction`       | Average number of mismatches on these reads per base                                           |
| `avg_sum_mismatch_qualities`           | Average sum of the base qualities of mismatches in the reads                                   |
| `num_q2_containing_reads`              | Number of reads with q2 runs at the 3' end                                                     |
| `avg_distance_to_q2_start_in_q2_reads` | Average distance of position (as fraction of unclipped read length) to the start of the q2 run |
| `avg_clipped_length`                   | Average clipped read length                                                                    |
| `avg_distance_to_effective_3p_end`     | Average distance to the 3' end of the read (as fraction of unclipped read length)              |
"""
function transformbamcounts(bamcounts::AbstractVector{String})
    # Declare an empty bam stats data frame
    countsdata = DataFrame(;
        chr                                  = String[],
        position                             = Int[],
        reference_base                       = String[],
        depth                                = Int[],
        base                                 = String[],
        count                                = Int[],
        avg_mapping_quality                  = Float64[],
        avg_basequality                      = Float64[],
        avg_se_mapping_quality               = Float64[],
        num_plus_strand                      = Int[],
        num_minus_strand                     = Int[],
        avg_pos_as_fraction                  = Float64[],
        avg_num_mismatches_as_fraction       = Float64[],
        avg_sum_mismatch_qualities           = Float64[],
        num_q2_containing_reads              = Int[],
        avg_distance_to_q2_start_in_q2_reads = Float64[],
        avg_clipped_length                   = Float64[],
        avg_distance_to_effective_3p_end     = Float64[],
    )

    # Transform the bam stats file
    for bamline in bamcounts
        # Split the base-independent stats by tabs
        bamfields = split(bamline, "\t")

        # Loop through the base-dependent stat blocks
        for i in 6:length(bamfields)
            # Split the base-dependent stats by colons
            basestats = split(bamfields[i], ":")

            # Parse the data into the correct types
            chr                                  = bamfields[1]
            position                             = parse(Int, bamfields[2])
            reference_base                       = bamfields[3]
            depth                                = parse(Int, bamfields[4])
            base                                 = basestats[1]
            count                                = parse(Int, basestats[2])
            avg_mapping_quality                  = parse(Float64, basestats[3])
            avg_basequality                      = parse(Float64, basestats[4])
            avg_se_mapping_quality               = parse(Float64, basestats[5])
            num_plus_strand                      = parse(Int, basestats[6])
            num_minus_strand                     = parse(Int, basestats[7])
            avg_pos_as_fraction                  = parse(Float64, basestats[8])
            avg_num_mismatches_as_fraction       = parse(Float64, basestats[9])
            avg_sum_mismatch_qualities           = parse(Float64, basestats[10])
            num_q2_containing_reads              = parse(Int, basestats[11])
            avg_distance_to_q2_start_in_q2_reads = parse(Float64, basestats[12])
            avg_clipped_length                   = parse(Float64, basestats[13])
            avg_distance_to_effective_3p_end     = parse(Float64, basestats[14])

            # Append the data to the dataframe
            push!(
                countsdata,
                [
                    chr,
                    position,
                    reference_base,
                    depth,
                    base,
                    count,
                    avg_mapping_quality,
                    avg_basequality,
                    avg_se_mapping_quality,
                    num_plus_strand,
                    num_minus_strand,
                    avg_pos_as_fraction,
                    avg_num_mismatches_as_fraction,
                    avg_sum_mismatch_qualities,
                    num_q2_containing_reads,
                    avg_distance_to_q2_start_in_q2_reads,
                    avg_clipped_length,
                    avg_distance_to_effective_3p_end,
                ],
            )
        end #for
    end #for

    return countsdata
end #function
