using HapLink

function clear_args()
    for i in 1:length(ARGS)
        pop!(ARGS)
    end #for
end #function

push!(ARGS, "variants")
push!(ARGS, "--bam")
push!(ARGS, "example/sample.bam")
push!(ARGS, "--reference")
push!(ARGS, "example/reference.fasta")
push!(ARGS, "--output")
push!(ARGS, "example/output.vcf")
push!(ARGS, "--quality")
push!(ARGS, "12")

HapLink.haplink()
clear_args()

push!(ARGS, "haplotypes")
push!(ARGS, "--bam")
push!(ARGS, "example/sample.bam")
push!(ARGS, "--variants")
push!(ARGS, "example/output.vcf")
push!(ARGS, "--output")
push!(ARGS, "example/output.yaml")
push!(ARGS, "--significance")
push!(ARGS, "0.05")
push!(ARGS, "--depth")
push!(ARGS, "50")
push!(ARGS, "--method")
push!(ARGS, "raw")

HapLink.haplink()

pop!(ARGS)
push!(ARGS, "ml-template")

HapLink.haplink()
clear_args()

push!(ARGS, "sequences")
push!(ARGS, "--haplotypes")
push!(ARGS, "example/output.yaml")
push!(ARGS, "--reference")
push!(ARGS, "example/reference.fasta")
push!(ARGS, "--output")
push!(ARGS, "example/output.fasta")

HapLink.haplink()
clear_args()

push!(ARGS, "--help")
HapLink.haplink()
