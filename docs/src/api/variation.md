# Variation Methods

```@meta
CurrentModule = HapLink
DocTestSetup = quote
    using HapLink
end
```

These methods extend the functionality of [`SequenceVariation.Variation`](https://github.com/BioJulia/SequenceVariation.jl)
for calculation of data related to `Variation`s created from NGS read alignments.

```@docs
seqpos
relativepos
quality(::Variation, ::Union{SAM.Record,BAM.Record})
variation(::VCF.Record, ::NucleotideSeq)
```
