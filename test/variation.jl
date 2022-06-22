const SUBSTITUTION = Variation(REFERENCE, "G8A")
const DELETION = Variation(REFERENCE, "Δ25-26")
const INSERTION = Variation(REFERENCE, "19GTC")

@testset "variation" begin
    @testset "variations" begin
        # Test that the types of variations is correct
        @test eltype(variations(VARIANTS)) <: Variation

        # Test for each type of Variation that should be present in the list
        @test SUBSTITUTION in variations(VARIANTS)
        @test DELETION in variations(VARIANTS)
        @test INSERTION in variations(VARIANTS)
    end #testset

    @testset "seqpos" begin
        # Test that the positions are what we expect
        @test seqpos(SUBSTITUTION, ALIGNMENTS[2]) == 8
        @test seqpos(DELETION, ALIGNMENTS[3]) == 24
        @test seqpos(INSERTION, ALIGNMENTS[end]) == 19
    end #testset

    @testset "relativepos" begin
        @test relativepos(SUBSTITUTION, SAMS[2]) == 8 / 33
        @test relativepos(DELETION, SAMS[3]) == 24 / 33
        @test relativepos(INSERTION, SAMS[end]) == 19 / 33
    end #testset

    @testset "quality" begin
        # The quality strings are constructed such that the quality corresponds to the
        # sequence position
        @test quality(SUBSTITUTION, SAMS[2]) == 8
        @test quality(DELETION, SAMS[3]) == 24.5
        @test quality(INSERTION, SAMS[end]) == 20
    end #testset
end #testset
