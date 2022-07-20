const VARIATIONPILEUP = VariationPileup(
    VARIATIONS[2], 0x04, [0.0, 1.0], [20.0, 30.0], [STRAND_POS, STRAND_NEG]
)

@testset "VariationPileup" begin
    @test altdepth(VARIATIONPILEUP) == 2
    @test frequency(VARIATIONPILEUP) == 0.5
    @test strand_bias(VARIATIONPILEUP) == 0.5
end #test
