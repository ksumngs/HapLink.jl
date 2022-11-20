@testset "VariationCall" begin
    @testset "call_variant" begin
        @test filters(call_variant(VARIATIONPILEUP, 1.0)) == ["PASS"]
    end #testset

    @testset "vcf_header" begin
        refpath = "/data/reference.fasta"
        α = 0.05
        D = 10
        Q = 12.0
        X = 0.1
        F = 0.05
        S = 0.3

        basebuffer = IOBuffer()
        write(basebuffer, HapLink._vcf_header(refpath, α))
        baseheader = String(take!(basebuffer))
        close(basebuffer)

        # Test the minimum header contents
        @test contains(baseheader, "##fileformat=VCFv4")
        @test contains(baseheader, "alpha=0.05")
        @test contains(baseheader, "##reference=file://$refpath")
        @test contains(baseheader, "##INFO=<ID=DP")

        # Test that filter fields _aren't_ There
        @test !contains(baseheader, "<ID=d")
        @test !contains(baseheader, "<ID=q")
        @test !contains(baseheader, "<ID=x")
        @test !contains(baseheader, "<ID=f")
        @test !contains(baseheader, "<ID=s")

        basebuffer = IOBuffer()
        write(basebuffer, HapLink._vcf_header(refpath, α; D=D, Q=Q, X=X, F=F, S=S))
        baseheader = String(take!(basebuffer))
        close(basebuffer)

        # Test that filter fields _are_ There
        @test contains(baseheader, "<ID=d")
        @test contains(baseheader, "<ID=q")
        @test contains(baseheader, "<ID=x")
        @test contains(baseheader, "<ID=f")
        @test contains(baseheader, "<ID=s")
    end #testset
end
