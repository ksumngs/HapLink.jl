<!-- markdownlint-disable -->

# <img src="./docs/src/assets/logo.png" style="border: 3px solid; float: left; margin: auto 2.5% auto 0" width="30%" > HapLink

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/ksumngs/HapLink.jl?label=version)](https://github.com/ksumngs/HapLink.jl/blob/master/CHANGELOG.md)
[![License](https://img.shields.io/github/license/ksumngs/HapLink.jl)](https://github.com/ksumngs/HapLink.jl/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ksumngs.github.io/HapLink.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ksumngs.github.io/HapLink.jl/dev)
[![Build Status](https://github.com/ksumngs/HapLink.jl/workflows/CI/badge.svg)](https://github.com/ksumngs/HapLink.jl/actions)
[![Coverage](https://codecov.io/gh/ksumngs/HapLink.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ksumngs/HapLink.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/haplink?color=green)](https://anaconda.org/bioconda/haplink)
[![Genie Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/HapLink)](https://pkgs.genieframework.com?packages=HapLink)

<!-- markdownlint-enable -->

> This project follows the [semver] _pro forma_ and uses the [OneFlow] branching
> model.

Call haplotypes based on linkage disequilibrium between variant sites on long
sequencing reads. Uses maximum liklihood methods for reads that are shorter than
the entire genome. Comes with its own variant caller.

---

## Installation

### :snake: Via Bioconda

**:warning::penguin: Linux-only!**

_Recommended for running HapLink on the **command line**_

```bash
conda create -n haplink -c bioconda -c conda-forge haplink -y
conda activate haplink
```

### ∴ Via Julia REPL

_Recommended for running HapLink within a **Julia session**_

```julia-repl
julia> ]
(@v1.6) pkg> add HapLink
```

To use this install of HapLink from the command line, you will need to add
`$HOME/.julia/bin` to your `$PATH`.

### :package: Via Apptainer

_Recommended for running HapLink on a **HPC**_

```bash
apptainer pull docker://ghcr.io/ksumngs/haplink.jl
```

### :whale: Via Docker

```bash
docker pull ghcr.io/ksumngs/haplink.jl:latest
```

## Usage

Please check [the docs] for more detailed instructions on how to use HapLink,
both on the command line and in the REPL.

The basic flow of HapLink is

<!-- markdownlint-disable -->

|     | Step                                   | Command              | Output Format |
| --- | -------------------------------------- | -------------------- | ------------- |
| 1.  | Call variants                          | `haplink variants`   | VCF           |
| 2.  | Call haplotypes                        | `haplink haplotypes` | YAML          |
| 3.  | Convert haplotype calls into sequences | `haplink sequences`  | FASTA         |

<!-- markdownlint-enable -->

You can see how this works using the files in the [example directory]:

```bash
haplink variants \
  example/reference.fasta \
  example/sample.bam \
  > sample.vcf

haplink haplotypes \
  example/reference.fasta \
  sample.vcf \
  example/sample.bam \
  > sample.yaml

haplink sequences \
  example/reference.fasta \
  sample.yaml \
  sample.fasta
```

## Development

HapLink is written in [Julia]. While the focus of the program is the command
line interface (CLI), it also exposes a nearly identical API in the form of a
Julia Package, which is described in [the docs].

### Development environment

For consistency, the recommended version of Julia as well as all the recommended
formatters and commit hooks are listed in a Nix file. If you have [direnv] and
[Nix] installed, then simply run

```bash
direnv allow .
pre-commit install
```

to setup Julia and the commit hook tools.

### Editing the package

HapLink.jl is a self-contained Julia package, and its development process is
identical to any other package as discussed in the [Pkg documentation].

```shellsession
$ git clone https://github.com/ksumngs/HapLink.jl.git
$ cd HapLink.jl
$ julia
(@v1.6) pkg> activate .
(HapLink) pkg> instantiate
julia> using HapLink
julia> ...
```

To test your changes on the command line application, ensure that
`$HOME/.julia/bin` is on your `$PATH`, then from the Julia REPL

```julia-repl
julia> ]
(@v1.6) pkg> activate .
(HapLink) pkg> build
```

This will update the application shim to include your changes.

[semver]: https://semver.org
[oneflow]:
  https://www.endoflineblog.com/oneflow-a-git-branching-model-and-workflow
[the docs]: https://ksumngs.github.io/HapLink.jl/stable
[example directory]: https://github.com/ksumngs/HapLink.jl/tree/master/example
[julia]: https://julialang.org
[pkg documentation]:
  https://pkgdocs.julialang.org/v1/managing-packages/#developing
