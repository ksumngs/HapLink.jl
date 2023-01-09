@testset "issam" begin
    @test HapLink._issam(SAM_FILE)
    @test HapLink._XAM_(SAM_FILE) == SAM
end

@testset "interval" begin
    @test interval(first(SAMS)) == Interval("REFERENCE", 1, 34)
end
