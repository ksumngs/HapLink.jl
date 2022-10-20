"""
    consensus(
        reference::Union{AbstractString,AbstractPath},
        variants::Union{AbstractString,AbstractPath};
        frequency::Float64=0.5,
        prefix::Union{AbstractString,Nothing}=nothing,
    )

Get the consensus `FASTA.Record `from `variants` applied to the first sequence in
`reference`.

# Arguments
- `reference::Union{AbstractString,AbstractPath}`: Path to the reference genome that
  variants were called from. Only the first sequence will be used.
- `variants::Union{AbstractString,AbstractPath}`: Path to variant call file that mutations
  will be applied from

# Keywords
- `frequency::Float64=0.5`: Fraction of total reads that must have supported the alternate
  position in order to be included as part of the consensus. In other words, only VCF
  records that have a `AF` (allele/alternate frequency) higher than this will be considered
  to contribute to the consensus.
- `prefix::Union{AbstractString,Nothing}=nothing`: Name to give to the output record. By
  default, the name of the output record will be the same as the name of the input record
  with `_CONSENSUS` appended. If `prefix` is supplied, then the name of the output record
  will be `\$(prefix)_CONSENSUS`.
"""
function consensus(
    reference::Union{AbstractString,AbstractPath},
    variants::Union{AbstractString,AbstractPath};
    frequency::Float64=0.5,
    prefix::Union{AbstractString,Nothing}=nothing,
)
    refrec = _first_record(reference)
    ref_id = isnothing(prefix) ? FASTA.identifier(refrec) : prefix
    ref_seq = FASTA.sequence(LongDNA{2}, refrec)

    con_seq = consensus(ref_seq, variants; frequency=frequency)

    fasta_record = FASTA.Record("$(ref_id)_CONSENSUS", con_seq)

    return fasta_record
end #function

"""
    consensus(
        reference::NucleotideSeq,
        variants::Union{AbstractString,AbstractPath};
        frequency::Float64=0.5,
    )

Get the consensus `FASTA.Record `from `variants` applied to the first sequence in
`reference`.

# Arguments
- `reference::NucleotideSeq`: Sequence of the reference genome that variants were called
  from
- `variants::Union{AbstractString,AbstractPath}`: Path to variant call file that mutations
  will be applied from

# Keywords
- `frequency::Float64=0.5`: Fraction of total reads that must have supported the alternate
  position in order to be included as part of the consensus. In other words, only VCF
  records that have a `AF` (allele/alternate frequency) higher than this will be considered
  to contribute to the consensus.
"""
function consensus(
    reference::NucleotideSeq,
    variants::Union{AbstractString,AbstractPath};
    frequency::Float64=0.5,
)
    SeqType = typeof(reference)
    BaseType = eltype(reference)

    vars = Variation{SeqType,BaseType}[]

    vcf_reader = VCF.Reader(open(string(variants), "r"))
    for vcf_rec in vcf_reader
        parse(Float64, VCF.info(vcf_rec, "AF")) >= frequency || continue
        all(f -> f == "PASS", VCF.filter(vcf_rec)) || continue
        push!(vars, variation(vcf_rec, reference))
    end #for

    con_seq = isempty(vars) ? reference : reconstruct!(reference, Variant(reference, vars))

    return con_seq
end #function
