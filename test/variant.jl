@test chromosome(var1) == "ref"
@test frequency(var1) == 6 / 10
@test frequency(var2) == 7 / 10

@test var1 ==
    Variant("ref", 10, "var1", dna"G", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
@test var1 !=
    Variant("fer", 10, "var1", dna"G", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
@test var1 !=
    Variant("ref", 11, "var1", dna"G", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
@test var1 !=
    Variant("ref", 10, "var2", dna"G", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
@test var1 !=
    Variant("ref", 10, "var1", dna"C", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
@test var1 !=
    Variant("ref", 10, "var1", dna"G", dna"T", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
@test var1 !=
    Variant("ref", 10, "var1", dna"G", dna"A", 29, :PASS, Dict(["DP" => 10, "AD" => 6]))
@test var1 !=
    Variant("ref", 10, "var1", dna"G", dna"A", 30, :PASS, Dict(["DP" => 11, "AD" => 6]))
