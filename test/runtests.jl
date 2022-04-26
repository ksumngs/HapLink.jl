using BioSequences
using Documenter
using HapLink
using Test

DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true)

@testset "HapLink.jl" begin
    include("sequences.jl")
    doctest(HapLink)
end
