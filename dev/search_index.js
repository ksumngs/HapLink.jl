var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API Reference","title":"API Reference","text":"CurrentModule = HapLink\nDocTestSetup = quote\n    using HapLink\nend","category":"page"},{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [HapLink]","category":"page"},{"location":"api/#HapLink.Haplotype","page":"API Reference","title":"HapLink.Haplotype","text":"Haplotype(mutations::AbstractVector{Variant})\n\nCreate an object to describe a group of mutations that appear together. A single Variant can be passed to create a haplotype of a single variant.\n\n\n\n\n\n","category":"type"},{"location":"api/#HapLink.Variant","page":"API Reference","title":"HapLink.Variant","text":"Variant\n\nDescribes a genomic mutation\n\nThe Variant type is based upon the Variant Call Format v4.2 specification, albeit with imperfect compliance.\n\nFields\n\nchromosome::String: An identifier from the reference genome or an angle-bracketed   ID String\nposition::Int: The reference position\nidentifier::String: Semicolon-separated list of unique identifiers where available. If   there is no identifier available, then \".\"\" value should be used.\nreferencebase::NucleotideSeq: Base at (or before, in case of insertion/deletion)   position in the reference genome\nalternatebase::NucleotideSeq: Base this mutation describes at position. Note that   each non-reference allele must be represented by a new Variant, unlike the VCF spec\nquality::Number: PHRED-scaled quality score for the assertion made by alternatebase\nfilter::Symbol: Filter status, :PASS is this position has passed all filters. Does Not   yet support multiple filters\ninfo::Dict{String,Any}: Additional information. No validation is made concerning the   keys or values.\n\nConstructors\n\nVariant(chromosome::String, position::Int, identifier::String, referencebase::NucleotideSeq,\n    alternatebase::NucleotideSeq, quality::Number, filter::Symbol, info::Dict{String,Any})\n\nVariant(data::DataFrameRow)\n\nVariants can be created from the default constructor, a VCF formatted string, or via a row generated by transformbamcounts.\n\nSee also Haplotype\n\n\n\n\n\n","category":"type"},{"location":"api/#HapLink.baseatreferenceposition-Tuple{XAM.BAM.Record, Int64}","page":"API Reference","title":"HapLink.baseatreferenceposition","text":"baseatreferenceposition(record::BAM.Record, pos::Int)\n\nGet the base at reference position pos present in the sequence of record. Returns nothing upon an error.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.callvariants-Tuple{DataFrames.AbstractDataFrame, Int64, Int64, Float64, Float64, Float64}","page":"API Reference","title":"HapLink.callvariants","text":"callvariants(bamcounts::AbstractDataFrame, D_min::Int, Q_min::Int, x_min::Float64,\n    f_min::Float64, α::Float64)\n\nBased on the aligned basecalls and stats in bamcounts, call variants.\n\nArguments\n\nbamcounts::AbstractDataFrame: DataFrame containing the output from bam-readcount\nD_min::Int: minimum variant depth\nQ_min::Int: minimum average PHRED-scaled quality at variant position\nx_min::Float64: minimum average fractional distance from read end at variant position\nf_min::Float64: minimum frequency of variant\nα::Float64: significance level of variants by Fisher's Exact Test\n\nReturns\n\nVector{Variant}: Variants that passed all the above filters\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.consensus-Tuple{BioSequences.LongDNASeq, AbstractArray{Variant}}","page":"API Reference","title":"HapLink.consensus","text":"consensus(refseq::LongDNASeq, vars::AbstractArray{Variant}; freq::Float64=0.5)\n\nGenerate the consensus sequence of refseq when mutated by vars, excluding any items from vars that have a variant frequency less than freq.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.containsposition-Tuple{XAM.BAM.Record, Int64}","page":"API Reference","title":"HapLink.containsposition","text":"containsposition(record::BAM.Record, i::Int)\n\nCheck to see if reference posision i is available in record\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.countbasestats-Tuple{String, String}","page":"API Reference","title":"HapLink.countbasestats","text":"countbasestats(bamfile::String, reffile::String)\n\nCount and calculate statistics on the basecalls of the alignment in bamfile to the reference genome in reffile. Returns a DataFrame with stats on every base in every alignment position. See transformbamcounts for a complete description of the output DataFrame schema.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.find_haplotypes-Tuple{AbstractVector{Variant}, AbstractString, AbstractString, Int64, Float64, Any}","page":"API Reference","title":"HapLink.find_haplotypes","text":"find_haplotypes(\n    variants::AbstractVector{Variant},\n    bamfile::AbstractString,\n    reffile::AbstractString,\n    D::Int,\n    α::Float64,\n    haplotypemethod\n)\n\nFind all combinations of variants that the reads in bamfile support as a valid haplotype with a minimum depth of D and Χ-squared linkage disequilibrium significance of α, with the haplotypes being converted into genomes via haplotypemethod.\n\nArguments\n\nvariants::AbstractVector{Variant}: A Vector of Variant objects which will be   combined into Haplotype objects and tested for validity\nbamfile::AbstractString: The path to a BAM file containing reads to check for the   presence of haplotypes\nreffile::AbstractString: The path to the reference genome in FASTA format\nD::Int: The minimum number of times a haplotype must be present in the reads according   to haplotypemethod\nα::Float64: The maximum Χ-squared p-value at which to consider a haplotype   significant and worth returning.\nhaplotypemethod: A function handle with the signature f(h::Haplotype,   b::AbstractString) which with return a table of the variant basecall matches found.   Both longread_genome and simulate_genome fulfil this requirement.\n\nReturns\n\nDict{Haplotype,Matrix{Int}}: A dictionary containing every significant haplotype and its   incidence matrix\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.linkage-Tuple{AbstractArray{Int64}}","page":"API Reference","title":"HapLink.linkage","text":"linkage(counts::AbstractArray{Int})\n\nCalculates the linkage disequilibrium and Chi-squared significance level of a combination of haplotypes whose number of occurrences are given by counts.\n\ncounts is an N-dimensional array where the Nth dimension represents the Nth variant call position within a haplotype. findoccurrences produces such an array.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.longread_genome-Tuple{Haplotype, AbstractString}","page":"API Reference","title":"HapLink.longread_genome","text":"longread_genome(haplotype::Haplotype, bamfile::AbstractString)\n\nParse the whole-genome length reads in bamfile to determine each read's basecall at every variant position within haplotype.\n\nArguments\n\nhaplotype::Haplotype: The combination of variants to test basecalls against\nbamfile::AbstractString: The path to a BAM file containing aligned reads to be tested   for instances of haplotype\n\nReturns\n\nMxN Array{Symbol} where M is the number of reads present in bamfile and   N=length(haplotype.mutations): A table of which base each read has in every variant   position of haplotype. The table has reads for rows and variant positions for columns,   e.g. longread_genome(...)[5,2] gives the basecall for the fifth read at the second   variant position. Basecalls are given as Symbol objects with possible values of\n:reference\n:alternate\n:other\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.matchvariant-Tuple{BioSequences.NucleotideSeq, Variant}","page":"API Reference","title":"HapLink.matchvariant","text":"matchvariant(base::Union{NucleotideSeq,DNA,AbstractVector{DNA}}, var::Variant)\n\nChecks if base matches the reference or variant expected in var, and returns a symbol indicating which, if any, it matches.\n\nReturned values can be :reference for a reference match, :alternate for an alternate match, or :other for no match with the given variant.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.mutate-Tuple{FASTX.FASTA.Record, Haplotype}","page":"API Reference","title":"HapLink.mutate","text":"function mutate(seq::FASTA.Record, haplotype::Haplotype)\nfunction mutate(seq::NucleotideSeq, haplotype::Haplotype)\n\nGive the reference sequence seq mutations to match the position and basecalls of haplotype. Returns a new sequence, leaving seq unmodified.\n\nWhen mutating a FASTA.Record, the new record is given a new unique identifier and description based on the SHA1 hash of the complete genotype.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.myref2seq-Tuple{BioAlignments.Alignment, Int64}","page":"API Reference","title":"HapLink.myref2seq","text":"myref2seq(aln::Alignment, i::Int)\n\nReplicates the functionality of BioAlignments ref2seq, but can handle hard clips by effectively removing them for the intent of finding the position.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.occurrence_matrix-Tuple{AbstractArray{Symbol}}","page":"API Reference","title":"HapLink.occurrence_matrix","text":"occurrence_matrix(readmatches::AbstractArray{Symbol})\n\nTransforms the haplotype occurrence table readmatches into an incidence matrix\n\nArguments\n\nreadmatches::AbstractArray{Symbol}: An mxn array where m is the number of   reads represented, and n is the number of variants in the haplotype considered, e.g.   readmatches[4,3] represents the match value for the third variant in the fourth read.   Valid values in the array are :reference, :alternate, and :other.\n\nReturns\n\n2x2x... Array{Int64, size(readmatches)[2]}: An N-dimensional matrix where N is   the number of variant positions in readmatches. The 1 index position in the   nth dimension represents the number of times the nth variant position was found   to have the reference base called, while the 2 index position represents the number   of times the nth variant position was found to have the alternate base called. E.g.   first(occurrence_matrix(reads)) gives the number of times the all-reference base   haplotype was found in reads, while occurrence_matrix(reads)[end] gives the number   of times the all-alternate base haplotype was found.\n\nExample\n\njulia> pseudoreads = [\n           :reference :reference :reference\n           :reference :reference :alternate\n           :reference :reference :alternate\n           :reference :reference :other\n       ];\n\njulia> occurrence_matrix(pseudoreads)\n2×2×2 Array{Int64, 3}:\n[:, :, 1] =\n 1  0\n 0  0\n\n[:, :, 2] =\n 2  0\n 0  0\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.overlap_inrange-Tuple{XAM.BAM.Record, XAM.BAM.Record}","page":"API Reference","title":"HapLink.overlap_inrange","text":"overlap_inrange(record1::BAM.Record, record2::BAM.Record; min::Int=0, max::Int=100)\n\nCheck if record1 and record2 overlap by an amount between min and max.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.phrederror-Tuple{Number}","page":"API Reference","title":"HapLink.phrederror","text":"phrederror(quality::Number)\n\nConverts a PHRED33-scaled error number into the expected fractional error of basecall\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.read_vcf-Tuple{Union{AbstractString, FilePathsBase.AbstractPath}}","page":"API Reference","title":"HapLink.read_vcf","text":"function read_vcf(vcf::AbstractPath)\n\nConvert the variant calls in vcf into Variants\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.savevcf-Tuple{AbstractVector{Variant}, AbstractString, AbstractString, Int64, Number, Float64, Float64}","page":"API Reference","title":"HapLink.savevcf","text":"savevcf(vars::AbstractVector{Variant}, savepath::String, refpath::String, D::Int,\n    Q::Number, x::Float64, α::Float64)\n\nSave a VCF file populated with vars\n\nArguments\n\nvars::AbstractVector{Variant}: Vector of Variants to write to file\nsavepath::AbstractString: path of the VCF file to write to. Will be overwritten\nrefpath::AbstractString: path of the reference genome used to call variants. The   absolute path will be added to the ##reference metadata\nD::Int: mimimum variant depth used to filter variants. Will be added as ##FILTER   metadata\nQ::Number: minimum PHRED quality used to filter variants. Will be added as ##FILTER   metadata\nx::Float64: minimum fractional read position used to filter variants. Will be added as   ##FILTER metadata\nα::Float64: Fisher's Exact Test significance level used to filter variants. Will be   added as ##FILTER metadata\n\nSaves the variants in vars to a VCF file at savepath, adding the reference genome refpath, the depth cutoff D, the quality cutoff Q, the position cutoff x, and the significance cutoff α as metadata.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.serialize_vcf-Tuple{Variant}","page":"API Reference","title":"HapLink.serialize_vcf","text":"serialize_vcf(v::Variant)\n\nCreate a VCF line to represent v.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.serialize_yaml-Tuple{Variant}","page":"API Reference","title":"HapLink.serialize_yaml","text":"serialize_yaml(v::Variant)\n\nCreate a valid YAML representation of v.\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.simulate_genome-Tuple{Haplotype, AbstractString}","page":"API Reference","title":"HapLink.simulate_genome","text":"simulate_genome(haplotype::Haplotype, bamfile::AbstractString; iterations::Int64=1000)\n\nSimulate a set of iterations reads that contain all of the variant positions in haplotype containing the actual reads present in bamfile via a maximum likelihood method.\n\nsimulate_genome examines each variant position in haplotype and finds a random read containing that position from bamfile. It will then check if the next variant position is contained on the previous read, and if not pick a new random read that contains that variant position. In this way, it assembles a set of reads that conceivably could have come from the same template strand via maximum likelihood.\n\nArguments\n\nhaplotype::Haplotype: The combination of variants to test the aligned reads for evidence   of\nbamfile::AbstractString: The path to a BAM file containing aligned reads to be tested   for instances of haplotype\n\nKeywords\n\niterations::Integer=1000: The number of times to combine reads and test for the presence   of haplotype\nnextreadcandidates: A function handle of the signature f(r1::BAM.Record,   r2::BAM.Record, pos::AbstractVecOrMat{Int}) that determines if r2 can be paired in   the same simulated genome as r2 based on the variant positions listed in pos\n\nReturns\n\nMxN Array{Symbol} where M=iterations and N=length(haplotype.mutations): A table of   which base each simulated read has in every variant position of haplotype. The table   has reads for rows and variant positions for columns, e.g. simulate_genome(...)[5,2]   gives the basecall for the fifth simulated read at the second variant position.   Basecalls are given as Symbol objects with possible values of\n:reference\n:alternate\n:other\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.sumsliced","page":"API Reference","title":"HapLink.sumsliced","text":"sumsliced(A::AbstractArray, dim::Int, pos::Int=1)\n\nSum all elements that are that can be referenced by pos in the dim dimension of A.\n\nExample\n\njulia> A = reshape(1:8, 2, 2, 2)\n2×2×2 reshape(::UnitRange{Int64}, 2, 2, 2) with eltype Int64:\n[:, :, 1] =\n 1  3\n 2  4\n\n[:, :, 2] =\n 5  7\n 6  8\n\njulia> sumsliced(A, 2)\n14\n\njulia> sumsliced(A, 2, 2)\n22\n\nHeavily inspired by Holy, Tim \"Multidimensional algorithms and iteration\" https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes\n\n\n\n\n\n","category":"function"},{"location":"api/#HapLink.transformbamcounts-Tuple{AbstractVector{String}}","page":"API Reference","title":"HapLink.transformbamcounts","text":"transformbamcounts(bamcounts::AbstractVector{String})\n\nConvert the output from bam-readcount to a DataFrame.\n\nSchema\n\nColumn name Description\nchr Chromosome\nposition Position\nreference_base Reference base\ndepth Total read depth at position\nbase Alternate base\ncount Number of reads containing base at position\navg_mapping_quality Mean mapping quality\navg_basequality Mean base quality at position\navg_se_mapping_quality Mean single-ended mapping quality\nnum_plus_strand Number of reads on the forward strand (N/A)\nnum_minus_strand Number of reads on the reverse strand (N/A)\navg_pos_as_fraction Average position on the read as a fraction, 0 (end) to 1 (center)\navg_num_mismatches_as_fraction Average number of mismatches on these reads per base\navg_sum_mismatch_qualities Average sum of the base qualities of mismatches in the reads\nnum_q2_containing_reads Number of reads with q2 runs at the 3' end\navg_distance_to_q2_start_in_q2_reads Average distance of position (as fraction of unclipped read length) to the start of the q2 run\navg_clipped_length Average clipped read length\navg_distance_to_effective_3p_end Average distance to the 3' end of the read (as fraction of unclipped read length)\n\n\n\n\n\n","category":"method"},{"location":"api/#HapLink.variant_positions_match-Tuple{XAM.BAM.Record, XAM.BAM.Record, AbstractVecOrMat{Int64}}","page":"API Reference","title":"HapLink.variant_positions_match","text":"variant_positions_match(\n    record1::BAM.Record,\n    record2::BAM.Record,\n    variantpositions::AbstractVecOrMat{Int}\n)\n\nCheck record1 and record2 at and return true if the basecalls are identical at every position in variantpositions.\n\n\n\n\n\n","category":"method"},{"location":"#HapLink","page":"Home","title":"HapLink","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for HapLink.","category":"page"}]
}
