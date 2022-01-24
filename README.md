
# <img src="./docs/src/assets/logo.png" style="border: 3px solid; float: left; margin: auto 2.5% auto 0" width="30%" > HapLink

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/ksumngs/HapLink.jl?label=version)](https://github.com/ksumngs/HapLink.jl/blob/master/CHANGELOG.md)
[![License](https://img.shields.io/github/license/ksumngs/HapLink.jl)](https://github.com/ksumngs/HapLink.jl/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ksumngs.github.io/HapLink.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ksumngs.github.io/HapLink.jl/dev)
[![Build Status](https://github.com/ksumngs/HapLink.jl/workflows/CI/badge.svg)](https://github.com/ksumngs/HapLink.jl/actions)
[![Coverage](https://codecov.io/gh/ksumngs/HapLink.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ksumngs/HapLink.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

> This project follows the [semver] _pro forma_ and uses the [OneFlow]
> branching model.

Call haplotypes based on linkage disequilibrium between variant sites on long
sequencing reads. Uses maximum liklihood methods for reads that are shorter than
the entire genome. Comes with its own variant caller.

***

## Installation

To run HapLink, your system must meet the following requirements

- Linux OS
- glibc
- x86-64 CPU

These restrictions apply even when using the package from within Julia.

If you need to use HapLink somewhere else, everything needed is available in a
[Docker image over on Quay].

### Prebuilt Binaries

To install our Hot-and-Ready binaries, run the following command:

```bash
mkdir -p ~/.local/opt/HapLink-0.5.1
curl -L https://github.com/ksumngs/HapLink.jl/releases/download/v0.5.1/HapLink-v0.5.1_linux.x86_64.tar.gz | tar xzv -C ~/.local/opt/HapLink-0.5.1
ln -s ~/.local/opt/HapLink-0.5.1/bin/haplink ~/.local/bin
```

### Julia Package

HapLink is not in the General Registry (yet!), so install using the `URL#tag`
syntax to use in the REPL.

```julia
using Pkg; Pkg.add("https://github.com/ksumngs/HapLink#v0.5.1")
```

## Usage

Please check [the docs] for actually useful instructions on how to use HapLink,
both on the command line and in the REPL.

The basic flow of HapLink is

|     | Step                                   | Command              | Output Format |
| --- | -------------------------------------- | -------------------- | ------------- |
| 1.  | Call variants                          | `haplink variants`   | VCF           |
| 2.  | Call haplotypes                        | `haplink haplotypes` | YAML          |
| 3.  | Convert haplotype calls into sequences | `haplink sequences`  | FASTA         |

You can see how this works using the files in the [example directory]:

```bash
haplink variants \
    --bam example/sample.bam \
    --reference example/reference.fasta
    --output sample.vcf

haplink haplotypes \
    --bam example/sample.bam \
    --variants sample.vcf \
    --output sample.yaml \

haplink sequences \
    --haplotypes sample.yaml \
    --reference example/reference.fasta \
    --output sample.fasta
```

## Development

HapLink is written in [Julia]. While the focus of the program is the command
line interface (CLI), it also exposes a nearly identical API in the form of a
Julia Package, which is described in [the docs].

### Editing the package

HapLink.jl is a self-contained Julia package, and its development process is
identical to any other package as discussed in the [Pkg documentation].
Personally, I tend to avoid the `dev` mode, and work straight from the cloned
package directory.

```shellsession
$ git clone https://github.com/ksumngs/HapLink.jl.git
$ cd HapLink.jl
$ julia
(@v1.6) pkg> activate .
(HapLink) pkg> instantiate
julia> using HapLink
julia> ...
```

### Creating the CLI application

To work with the CLI directly, you can do one (or both) of the following

#### 1. Create a shim

> - _Fast to implement_
> - _Changes are reflected immediately_
> - _Slow execution time ([TTFP])_

In my `~/bin` directory, I have an executable file named `haplink` with the
following contents:

```bash
#!/bin/sh
julia --project=$HOME/src/HapLink.jl -e 'using HapLink.haplink()' "$@"
```

#### 2. Compile the binary

> - _More involved implementation_
> - _Updates must be recompiled_
> - _Fast execution time_

Binaries are compiled using [PackageCompiler.jl], using the recipe in [.github/workflows/build.yml].

1. Get the [official Julia release] (disto packages generally don't work)
2. Install PackageCompiler into that Julia depot

    ```shellsession
    (@v1.6) pkg> install PackageCompiler
    ```

3. Run `PackageCompiler.create_app()` with the following options

    ```julia
    using PackageCompiler
    create_app(
      "/path/to/HapLink.jl",
      "/path/to/output",
      precompile_execution_file="precompile_app.jl",
      executables=["haplink" => "haplink"],
      cpu_target="x86-64",
    )
    ```

Compilation can take over 15 minutes to complete, so be patient!

## Contributors

It's pretty lonely here: HapLink was solely made by Thomas Christensen while
working at Kansas State University. Why don't you [open a pull request] and fix
that?

[semver]: https://semver.org
[OneFlow]: https://www.endoflineblog.com/oneflow-a-git-branching-model-and-workflow
[Docker image over on Quay]: https://quay.io/repository/millironx/julia_bam-readcounts
[the docs]: https://ksumngs.github.io/HapLink.jl/stable
[example directory]: https://github.com/ksumngs/HapLink.jl/tree/master/example
[Julia]: https://julialang.org
[Pkg documentation]: https://pkgdocs.julialang.org/v1/managing-packages/#developing
[TTFP]: https://viralinstruction.com/posts/badjulia/#compile_time_latency
[PackageCompiler.jl]: https://julialang.github.io/PackageCompiler.jl/stable/apps.html
[.github/workflows/build.yml]: https://github.com/ksumngs/HapLink.jl/blob/master/.github/workflows/build.yml
[official Julia release]: https://julialang.org/downloads/
[open a pull request]: https://github.com/ksumngs/HapLink.jl/compare
