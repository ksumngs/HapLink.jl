
# <img src="./docs/src/assets/logo.png" style="border: 3px solid; float: left; margin: auto 2.5% auto 0" width="30%" > HapLink

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/ksumngs/HapLink.jl?label=version)](https://github.com/ksumngs/HapLink.jl/blob/master/CHANGELOG.md)
[![License](https://img.shields.io/github/license/ksumngs/HapLink.jl)](https://github.com/ksumngs/HapLink.jl/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ksumngs.github.io/HapLink.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ksumngs.github.io/HapLink.jl/dev)
[![Build Status](https://github.com/ksumngs/HapLink.jl/workflows/CI/badge.svg)](https://github.com/ksumngs/HapLink.jl/actions)
[![Coverage](https://codecov.io/gh/ksumngs/HapLink.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ksumngs/HapLink.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

> This project follows the [semver] _pro forma_ and uses the [git-flow]
> branching model.

Call haplotypes based on linkage disequilibrium between variant sites on long
sequencing reads. Uses maximum liklihood methods for reads that are shorter than
the entire genome. Comes with its own variant caller.

***

## Installation

You must have both [Julia] version 1.6 or higher and [bam-readcount] version 1.0
or higher installed and available on your `PATH`. Put down your pitchforks:
everything needed is available in a [Docker image over on Quay] until the
bam-readcount dependency can be removed.

To use as a command-line program on Unix systems, download the [latest release]
to a safe place, then symlink the `haplink` file to somewhere on your `PATH`.
This method is not available on Windows.

```bash
wget -qO- https://github.com/ksumngs/HapLink.jl/archive/refs/tags/v0.1.0.tar.gz | tar xvz -C $HOME/.local/opt
ln -s $HOME/.local/opt/HapLink.jl-0.1.0/haplink $HOME/.local/bin
```

HapLink is not in the General Registry (yet!), so install using the `URL#tag`
syntax to use in the REPL. This works on all operating systems that Julia does.

```julia
using Pkg; Pkg.add("https://github.com/ksumngs/HapLink#v0.1.0")
```

## Usage

Please check [the docs] for actually useful instructions on how to use HapLink,
both on the command line and in the REPL.

You will need to provide a BAM file containing aligned reads and a FASTA file
containing the reference genome those reads were aligned to. (Note that the
reference **is** mandatory, while the help text might indicate otherwise. See
https://github.com/carlobaldassi/ArgParse.jl/issues/108.)

HapLink will output

1. A YAML file containing metadata on all found haplotypes
2. A FASTA file containing the sequences of all found haplotypes
3. (Optional) A VCF file containing the variant calls used to call haplotypes

```shellsession
$ haplink --help
usage: haplink -r REFERENCE [-g ANNOTATIONS] [-v VARIANTS] [-p PREFIX]
               [-q QUALITY] [-t FREQUENCY] [-x POSITION]
               [-a VARIANT_SIGNIFICANCE] [-b HAPLOTYPE_SIGNIFICANCE]
               [-d VARIANT_DEPTH] [-u HAPLOTYPE_DEPTH]
               [--method METHOD] [--iterations ITERATIONS] [--version]
               [-h] bamfile

A haplotype caller for long sequencing reads using linkage
disequilibrium

positional arguments:
  bamfile               BAM formatted file containing the alignment to
                        call haplotypes from

optional arguments:
  -r, --reference REFERENCE
                        FASTA formatted reference genome
  -g, --annotations ANNOTATIONS
                        GFF3 formatted annotations for the reference
                        genome (NOT YET IMPLEMENTED)
  -v, --variants VARIANTS
                        File to output all variant calls in VCF
                        format
  -p, --prefix PREFIX   Test to start the output file names with. If
                        unspecified, will use the name of the
                        alignment file up to the first dot (.)
  -q, --quality QUALITY
                        Minimum average quality (PHRED score) of a
                        variant basecall (type: Int64, default: 21)
  -t, --frequency FREQUENCY
                        Minimum frequency at which a variant must
                        appear (type: Float64, default: 0.05)
  -x, --position POSITION
                        Remove variants that occur only in positions
                        within this percentage of the end (type:
                        Float64, default: 0.1)
  -a, --variant-significance VARIANT_SIGNIFICANCE
                        Maximum Χ-squared significance level to
                        consider a haplotype (type: Float64, default:
                        1.0e-5)
  -b, --haplotype-significance HAPLOTYPE_SIGNIFICANCE
                        Maximum Χ-squared significance level to
                        consider a haplotype (type: Float64, default:
                        1.0e-5)
  -d, --variant-depth VARIANT_DEPTH
                        Minimum depth to consider a variant (type:
                        Int64, default: 10)
  -u, --haplotype-depth HAPLOTYPE_DEPTH
                        Minimum depth to consider a haplotype (type:
                        Int64, default: 10)
  --method METHOD       Haplotype read-building method. Choose one of
                        'ml-overlap', 'ml-gapped' or 'raw' (NOT YET
                        IMPLEMENTED) (default: "ml-overlap")
  --iterations ITERATIONS
                        Formula to determine how many iterations to
                        perform using one of the ml methods (NOT YET
                        IMPLEMENTED) (default: ":(1000)")
  --version             show version information and exit
  -h, --help            show this help message and exit
```

### Example

```shellsession
$ haplink /data/sample.bam -r /data/reference.fasta -v variants.vcf
Minimum mapping quality is set to 0
$ ls
sample.fasta
sample.vcf
sample.yaml
```

## Contributors

It's pretty lonely here: HapLink was solely made by Thomas Christensen while
working at Kansas State University. Why don't you [open a pull request] and fix
that?

[bam-readcount]: https://github.com/genome/bam-readcount
[Docker image over on Quay]: https://quay.io/repository/millironx/julia_bam-readcounts
[git-flow]: https://nvie.com/posts/a-successful-git-branching-model
[Julia]: https://julialang.org
[latest release]: https://github.com/ksumngs/HapLink.jl/releases/latest
[open a pull request]: https://github.com/ksumngs/HapLink.jl/compare
[semver]: https://semver.org
[the docs]: https://ksumngs.github.io/HapLink.jl/stable
