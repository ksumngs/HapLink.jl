@testset "VariationCall" begin
    @test filters(call_variant(VARIATIONPILEUP, 1.0)) == ["PASS"]
end
