using HapLink
using Test
using Documenter

DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true)

@testset "HapLink.jl" begin
    doctest(HapLink)
end
