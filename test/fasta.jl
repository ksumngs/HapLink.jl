@testset "fasta" begin
    fasta_file, _ = mktemp()
    open(fasta_file, "w") do f
        write(
            f,
            """>REFERENCE
$REFERENCE
>GENOTYPE_V
$GENOTYPE_V
   """,
        )
    end #do
    @test FASTA.identifier(HapLink._first_record(fasta_file)) == "REFERENCE"
    @test FASTA.sequence(LongDNA{4}, HapLink._first_record(fasta_file)) == REFERENCE
end #testset
