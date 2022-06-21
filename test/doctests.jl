DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true)

@testset "doctests" begin
    doctest(HapLink)
end #testset
