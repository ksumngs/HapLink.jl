@testitem "Haplotype CIGAR" begin
    using BioSymbols: DNA
    using BioSequences: LongDNA, @dna_str
    using SequenceVariation: Haplotype, Variation

    reference = dna"TGATGCGTGTAGCAACACTTATAGCG"
    reference_genotype = Haplotype(reference, Variation{LongDNA{4},DNA}[])
    genotype = Haplotype(
        reference,
        [
            Variation(reference, "Δ1-2"),
            Variation(reference, "10T"),
            Variation(reference, "Δ17-18"),
            Variation(reference, "A23C"),
        ],
    )

    @test cigar(reference_genotype) == "26M"
    @test cigar(genotype) == "2D8M1I6M2D8M"
end #testitem
