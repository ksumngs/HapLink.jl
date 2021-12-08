using HapLink

push!(ARGS, "example/sample.bam")
push!(ARGS, "--reference")
push!(ARGS, "example/reference.fasta")
push!(ARGS, "--variants")
push!(ARGS, "sample.vcf")
push!(ARGS, "--quality")
push!(ARGS, "12")
push!(ARGS, "--haplotype-significance")
push!(ARGS, "0.05")
push!(ARGS, "--haplotype-depth")
push!(ARGS, "50")

HapLink.haplink()

push!(ARGS, "--method")
push!(ARGS, "raw")

HapLink.haplink()

for i in 1:length(ARGS)
    pop!(ARGS)
end #for

push!(ARGS, "example/sample.yaml")
push!(ARGS, "example/reference.fasta")
push!(ARGS, "haplotypes.fasta")

HapLink.make_haplotype_fastas()

for i in 1:length(ARGS)
    pop!(ARGS)
end #for

push!(ARGS, "--help")
HapLink.haplink()
