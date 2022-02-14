using BioSymbols
using GenomicFeatures
using VariantCallFormat
using XAM

@testset "SNPs" begin
    # Set up a correct SNP
    snp = SNP(Interval("X", 1, 1), DNA_A, DNA_T)

    # Check the accessor functions
    @test location(snp) == Interval("X", 1, 1)
    @test refbase(snp) == DNA_A
    @test altbase(snp) == DNA_T

    # Check that we don't mix DNA and RNA
    @test_throws MethodError SNP(Interval("X", 1, 1), DNA_T, RNA_U)

    # Check that we can create a snp from a VCF record
    @test SNP(VCF.Record("X\t1\t.\tA\tT\t30.0\tPASS\t")) == snp

    # Check that we can't create a snp from a multi-insertion VCF record
    @test_throws ArgumentError SNP(VCF.Record("X\t1\t.\tA\tTT\t30.0\tPASS\t"))

    # Check that RNA is properly parsed
    @test typeof(SNP(VCF.Record("X\t1\t.\tA\tU\t30.0\tPASS\t"))) == SNP{RNA}
end
