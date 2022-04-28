@test Haplotype([var1]) == Haplotype([
    Variant("ref", 10, "var1", dna"G", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
])
@test Haplotype([var1, var3]) == Haplotype([var3, var1])
@test Haplotype([var1, var3]) != Haplotype([var1])
