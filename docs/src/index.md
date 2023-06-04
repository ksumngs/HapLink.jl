# HapLink

```@meta
CurrentModule = HapLink
DocTestSetup = quote
    using HapLink
end
```

```@docs
HapLink
```

## Welcome

Howdy! 🤠 And welcome to HapLink! 👋 HapLink is a command-line suite of tools to
enable the exploration of viral quasispecies within a single metagenomic sample.
Every piece eventually builds up to our viral haplotype caller, which uses
linkage disequilibrium on long sequencing reads (💡 think
[Oxford Nanopore](https://nanoporetech.com/) or
[PacBio HiFi](https://www.pacb.com/)) to identify genetic mutations that are
conserved within a single virus particle.

This manual will cover the different ways of using HapLink, starting with a few
tutorials before diving into the details of our reference section.

### Contents

```@contents

```

## Getting started

Ready to dive in? 🤿 Here's a 30,000-foot view

```bash
curl -fsSL https://install.julialang.org | sh -s -- --yes
source ~/.bashrc
julia \
  --startup-file=no \
  --history-file=no \
  --quiet \
  -e 'using Pkg; Pkg.activate(;temp=true); Pkg.add(HapLink)'
echo 'export PATH=$HOME/.julia/bin:$PATH' >> $HOME/.bashrc
source ~/.bashrc

wget https://github.com/ksumngs/HapLink.jl/raw/v1.0.0-rc1/example/reference.fasta
wget https://github.com/ksumngs/HapLink.jl/raw/v1.0.0-rc1/example/sample.bam
wget https://github.com/ksumngs/HapLink.jl/raw/v1.0.0-rc1/example/sample.bam.bai

haplink variants \
  reference.fasta \
  sample.bam \
  --significance 0.5 \
  --depth 1 \
  --quality 10.0 \
  --position 0.01 \
  --frequency 0.05 \
  | tee sample.vcf

haplink consensus \
  reference.fasta \
  sample.vcf \
  | tee consensus.fasta

haplink haplotypes \
  reference.fasta \
  sample.vcf \
  sample.bam \
  --frequency 0.75 \
  | tee sample.yml

haplink sequences \
  reference.fasta \
  sample.yml \
  | tee sample.fasta
```
