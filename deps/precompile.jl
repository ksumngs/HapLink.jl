using HapLink

const PACKAGE_DIR = dirname(dirname(pathof(HapLink)))
const EXAMPLE_DIR = joinpath(PACKAGE_DIR, "example")

HapLink.variants(
    joinpath(EXAMPLE_DIR, "reference.fasta"),
    joinpath(EXAMPLE_DIR, "sample.bam");
    outfile=joinpath(EXAMPLE_DIR, "sample.vcf"),
    significance=0.5,
    depth=UInt(1),
    quality=10.0,
    position=0.01,
    frequency=0.05,
)

HapLink.consensus(
    joinpath(EXAMPLE_DIR, "reference.fasta"),
    joinpath(EXAMPLE_DIR, "sample.vcf");
    outfile=joinpath(EXAMPLE_DIR, "consensus.fasta"),
)

HapLink.haplotypes(
    joinpath(EXAMPLE_DIR, "reference.fasta"),
    joinpath(EXAMPLE_DIR, "sample.vcf"),
    joinpath(EXAMPLE_DIR, "sample.bam");
    outfile=joinpath(EXAMPLE_DIR, "sample.yml"),
    frequency=0.75,
)

HapLink.haplotypes(
    joinpath(EXAMPLE_DIR, "reference.fasta"),
    joinpath(EXAMPLE_DIR, "sample.vcf"),
    joinpath(EXAMPLE_DIR, "sample.bam");
    outfile=joinpath(EXAMPLE_DIR, "sample.ml.yml"),
    frequency=0.75,
    simulated_reads=true,
    iterations=UInt(10),
)

HapLink.sequences(
    joinpath(EXAMPLE_DIR, "reference.fasta"),
    joinpath(EXAMPLE_DIR, "sample.yml");
    outfile=joinpath(EXAMPLE_DIR, "sample.fasta"),
)
