const VARIATIONINFO = VariationInfo(VARIATIONS[1], 0.5, 30.0, STRAND_POS)

@testset "variationinfo" begin
    @testset "getters" begin
        @test variation(VARIATIONINFO) == VARIATIONS[1]
        @test readpos(VARIATIONINFO) == 0.5
        @test quality(VARIATIONINFO) == 30
        @test strand(VARIATIONINFO) == STRAND_POS
    end #testset

    @testset "variationinfos" begin
        @test isempty(variationinfos(SAMS[1], REFERENCE))
        @test SUBSTITUTION in variation.(variationinfos(SAMS[2], REFERENCE))
        @test DELETION in variation.(variationinfos(SAMS[3], REFERENCE))
        @test INSERTION in variation.(variationinfos(SAMS[end], REFERENCE))
    end #testset
end #testset
