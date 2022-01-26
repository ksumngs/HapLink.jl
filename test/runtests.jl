using BioSequences
using Documenter
using HapLink
using Test

const REFERENCE = dna"AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT"
var1 = Variant("ref", 10, "var1", dna"G", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
var2 = Variant("ref", 10, "var2", dna"G", dna"T", 30, :PASS, Dict(["DP" => 10, "AD" => 7]))
var3 = Variant("ref", 20, "var3", dna"C", dna"G", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))

DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true)

@testset "HapLink.jl" begin
    include("variant.jl")
    include("sequences.jl")
    include("haplotype.jl")
    doctest(HapLink)
    include("snptests.jl")
end
