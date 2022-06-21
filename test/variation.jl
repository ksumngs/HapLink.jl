@testset "variation" begin
    @testset "variations" begin
        # Test that the types of variations is correct
        @test eltype(variations(VARIANTS)) <: Variation

        # Test for each type of Variation that should be present in the list
        @test Variation(REFERENCE, "G8A") in variations(VARIANTS)
        @test Variation(REFERENCE, "Δ25-26") in variations(VARIANTS)
        @test Variation(REFERENCE, "19GUC") in variations(VARIANTS)
    end #testset

    @testset "seqpos" begin
        # Test that the positions are what we expect
        @test seqpos(VARIATIONS[1], ALIGNMENTS[2]) == 8
        @test seqpos(VARIATIONS[5], ALIGNMENTS[3]) == 24
        @test seqpos(VARIATIONS[end - 2], ALIGNMENTS[end]) == 19
    end #testset

    @testset "relativepos" begin
        @test relativepos(VARIATIONS[1], SAMS[2]) == 8 / 33
        @test relativepos(VARIATIONS[5], SAMS[3]) == 24 / 33
        @test relativepos(VARIATIONS[end - 2], SAMS[end]) == 19 / 33
    end #testset
end #testset
