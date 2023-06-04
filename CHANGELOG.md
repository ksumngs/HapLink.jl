# Changelog

All notable changes to HapLink.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

> This release represents a [big rewrite](https://youtu.be/xCGu5Z_vaps) of
> HapLink

### Added

- BAM index file finder ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- `VariationInfo` type ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- `VariationPileup` type ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- `VariationCall` type ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- `HaplotypeCall` type ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- More rigorous haplotype combination finding algorithm
  ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- `significance` function for haplotype chi-squared significance
  ([#37](https://github.com/ksumngs/HapLink.jl/pull/37))
- Earthfile for repeatable testing and documentation builds
  ([#42](https://github.com/ksumngs/HapLink.jl/issues/42)/[#45](https://github.com/ksumngs/HapLink.jl/pull/45))
- Nix configuration
  ([#41](https://github.com/ksumngs/HapLink.jl/issues/41)/[#53](https://github.com/ksumngs/HapLink.jl/pull/52))

### Changed

- VCF parsing outsourced to
  [VariantCallFormat.jl](https://github.com/rasmushenningsson/VariantCallFormat.jl)
  ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- Codebase now uses [Blue Style](https://github.com/invenia/BlueStyle)
  ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- Variant calling and storage outsourced to
  [SequenceVariation.jl](https://github.com/BioJulia/SequenceVariation.jl)
  ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- [BREAKING] `linkage` now returns only the unweighted linkage disequilibrium,
  not a tuple of linkage disequilibrium and significance
  ([#37](https://github.com/ksumngs/HapLink.jl/pull/37))
- CLI now implemented in [Comonicon.jl](https://comonicon.org) instead of
  [ArgParse](https://github.com/carlobaldassi/ArgParse.jl/)
  ([#38](https://github.com/ksumngs/HapLink.jl/issues/38)/[#47](https://github.com/ksumngs/HapLink.jl/pull/47))
- `haplink haplotypes` now only outputs simulation-related settings if
  `--simulated-reads` flag is set
  ([#40](https://github.com/ksumngs/HapLink.jl/issues/40)/[#48](https://github.com/ksumngs/HapLink.jl/pull/48))

### Removed

- [BREAKING] `myref2seq` (functionality now present in
  [BioAlignments.jl](https://github.com/BioJulia/BioAlignments.jl))
  ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- Dependency on bam-readcount artifact (including x86-64 glibc Linux-specific
  code) ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- `Variant` type (functionality now present in
  [SequenceVariation.jl](https://github.com/BioJulia/SequenceVariation.jl))
  ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))
- `Haplotype` type (functionality now present in
  [SequenceVariation.jl](https://github.com/BioJulia/SequenceVariation.jl))
  ([#35](https://github.com/ksumngs/HapLink.jl/pull/35))

### Fixed

- Linkage disequilibrium and chi-squared significance calculations corrected
  ([#37](https://github.com/ksumngs/HapLink.jl/pull/37))

## [0.7.1] - 2022-07-29

### Fixed

- VCF output files no longer contain invalid characters
  ([#30](https://github.com/ksumngs/HapLink.jl/pull/30))

## [0.7.0] - 2022-04-28

### Added

- `haplink consensus` command to generate consensus sequences from variant calls
  ([#29](https://github.com/ksumngs/HapLink.jl/pull/29))

### Changed

- `haplink haplotypes` now calls haplotypes based on the consensus sequence
  ([#30](https://github.com/ksumngs/HapLink.jl/pull/30))

## [0.6.1] - 2022-03-28

### Fixed

- The CI scripts have been fixed

## [0.6.0] - 2022-03-28

### Changed

- `haplink haplotypes` now requires a reference genome passed via `--reference`
  ([PR#25](https://github.com/ksumngs/HapLink.jl/pull/25))

## [0.5.1] - 2022-01-24

### Added

- `--seed` parameter for fixing random seed

## [0.5.0] - 2022-01-12

### Removed

- `bam-readcount` dependency

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

- New Julia Package via
  [PkgTemplates.jl](https://github.com/invenia/PkgTemplates.jl)
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

[unreleased]: https://github.com/ksumngs/HapLink.jl/compare/v0.7.1...HEAD
[0.7.1]: https://github.com/ksumngs/HapLink.jl/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/ksumngs/HapLink.jl/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.4.3...v0.5.0
[0.4.3]: https://github.com/ksumngs/HapLink.jl/compare/v0.4.2...v0.4.3
[0.4.2]: https://github.com/ksumngs/HapLink.jl/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/ksumngs/HapLink.jl/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/ksumngs/HapLink.jl/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/ksumngs/HapLink.jl/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/ksumngs/HapLink.jl/compare/v0.0.1...v0.1.1
[0.0.1]: https://github.com/ksumngs/yavsap/releases/tag/v0.0.1
