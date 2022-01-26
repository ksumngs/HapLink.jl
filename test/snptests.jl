using BioSymbols
using GenomicFeatures
using VariantCallFormat
using XAM

@testset "SNPs" begin
    # Set up a correct SNP
    snp = SNP(Interval("X", 1, 1), DNA_A, DNA_T)

    # Check that we don't mix DNA and RNA
    @test_throws MethodError SNP(Interval("X", 1, 1), DNA_T, RNA_U)

    # Check that we can create a snp from a VCF record
    @test SNP(VCF.Record("X\t1\t.\tA\tT\t30.0\tPASS\t")) == snp

    # Check that we can't create a snp from a multi-insertion VCF record
    @test_throws ArgumentError SNP(VCF.Record("X\t1\t.\tA\tTT\t30.0\tPASS\t"))

    # Check that RNA is properly parsed
    @test typeof(SNP(VCF.Record("X\t1\t.\tA\tU\t30.0\tPASS\t"))) == SNP{RNA}

    # Check that removing padded reads doesn't totally invalidate results
    samrecords = SAM.Record.([
        "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGACACTG\t*",
        "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t*",
        "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t*\tSA:Z:ref,29,-,6H5M,17,0;",
        "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tACAGCTTCAGC\t*",
        "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t*\tSA:Z:ref,9,+,5S6M,30,1;",
        "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t*\tNM:i:1"
    ])
    @test (@test_logs (:warn, """
    One or more of these reads contains a 'P' CIGAR operation, which is not yet
    supported by BioAlignments.jl. These reads will be ignored.
    """) depth(SNP("ref", 17, DNA_T, DNA_C), samrecords)) == 2

end
