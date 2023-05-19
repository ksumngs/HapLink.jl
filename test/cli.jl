@testset "CLI" begin
    @test include(joinpath(dirname(dirname(pathof(HapLink))), "deps", "precompile.jl")) == 0
end #testset
