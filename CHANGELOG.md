# Changelog

All notable changes to HapLink.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.3] - 2021-12-30

### Fixed

- Can now create and manipulate empty (reference) haplotypes

## [0.4.2] - 2021-12-29

### Changed

- Use CHANGELOG for release notes on GitHub

## [0.4.1] - 2021-12-28

### Changed

- Dockerfile now separated into multiple layers

### Fixed

- Precompile script updated to use HapLink 0.4 output format

## [0.4.0] - 2021-12-28

### Changed

- Output YAML schema changed to include more information

## [0.3.0] - 2021-12-28

### Changed

- HapLink CLI functionality is now exposed as subcommands:
  - `haplink variants`
  - `haplink haplotypes`
  - `haplink sequences`
- Dockerfile version numbers extracted to environment variables

### Removed

- `make-haplotype-fastas` command line/binary application
- `haplink haplotypes` no longer outputs in fasta format

## [0.2.0] - 2021-12-08

### Added

- Julia Formatter to CI workflow
- Example files
- Precompile script
- Long reads/raw read haplotype finding

### Changed

- Extacted the incidence matrix calculations to `occurrence_matrix` function
- Improved REPL printout for `Haplotype` type
- Haplotype read finder must now be passed as argument to haplotype finder

### Removed

- Mac and Windows binary builds

## [0.1.2] - 2021-12-07

### Changed

- Correct build workflow typo

## [0.1.1] - 2021-12-07

### Changed

- Moved sequence manipulation functions to `sequences.jl`
- Moved haplotype finding functions to `haplotypecalling.jl`
- Moved variant calling functions to `variantcalling.jl`
- Improved binary building workflow

## [0.0.1] - 2021-12-06

### Added

- New Julia Package via [PkgTemplates.jl](https://github.com/invenia/PkgTemplates.jl)
- GitHub Issue Templates
- Code of Conduct
- Self-contained Dockerfile
- Logo artwork
- CLI programs
  - `haplink`
  - `make-haplotype-fastas`
- Types
  - `Haplotype`
  - `Variant`

[0.4.3]: https://github.com/ksumngs/HapLink.jl/compare/v0.4.2...v0.4.3
[0.4.2]: https://github.com/ksumngs/HapLink.jl/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/ksumngs/HapLink.jl/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/ksumngs/HapLink.jl/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/ksumngs/HapLink.jl/compare/v0.0.1...v0.1.1
[0.0.1]: https://github.com/ksumngs/yavsap/releases/tag/v0.0.1