@testitem "Pseudoread lendiff" begin
    using HapLink: Pseudoread, _lendiff
    using SequenceVariation: Haplotype, Variation, mutation
    using BioSequences: LongDNA, @dna_str
    using BioSymbols: DNA

    const REFSEQ = dna"GATTACA"
    const SUB = Variation(REFSEQ, "T3G")
    const INS = Variation(REFSEQ, "4GATTACA")
    const DEL = Variation(REFSEQ, "Δ5-6")

    # Test that empty haplotype gives the same number of bases
    @test _lendiff(Pseudoread(1, 7, Haplotype(REFSEQ, Variation{LongDNA{4},DNA}[]))) == 0

    # Test that a haplotype with only substitutions gives the same number of bases
    @test _lendiff(Pseudoread(1, 7, Haplotype(REFSEQ, [SUB]))) == 0

    # Test that a haplotype with an insertion gives the correct number of extra bases
    @test _lendiff(Pseudoread(1, 7, Haplotype(REFSEQ, [INS]))) == 7

    # Test that a haplotype with a deletion gives the correct number of deleted bases
    @test _lendiff(Pseudoread(1, 7, Haplotype(REFSEQ, [DEL]))) == -2

    # Test that a combined haplotype gives the combined number of bases
    @test _lendiff(Pseudoread(1, 7, Haplotype(REFSEQ, [SUB, INS, DEL]))) == 5

    # Test that removing the insertion from the read window removes the added bases
    @test _lendiff(Pseudoread(1, 4, Haplotype(REFSEQ, [INS]))) == 0
    @test _lendiff(Pseudoread(1, 5, Haplotype(REFSEQ, [INS]))) == 7

    # Test that removing the deletion from the read window removes the deleted bases
    @test _lendiff(Pseudoread(1, 4, Haplotype(REFSEQ, [DEL]))) == 0

    # Test that _partially_ removing the deletion _lessens_ the number of deleted bases
    @test _lendiff(Pseudoread(1, 5, Haplotype(REFSEQ, [DEL]))) == -1
    @test _lendiff(Pseudoread(6, 7, Haplotype(REFSEQ, [DEL]))) == -1
end #@testitem
