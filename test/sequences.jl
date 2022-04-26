const REFERENCE = dna"AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT"
var1 = Variant("ref", 10, "var1", dna"G", dna"A", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))
var2 = Variant("ref", 10, "var2", dna"G", dna"T", 30, :PASS, Dict(["DP" => 10, "AD" => 7]))
var3 = Variant("ref", 20, "var3", dna"C", dna"G", 30, :PASS, Dict(["DP" => 10, "AD" => 6]))

# Test single point mutation
@test consensus(REFERENCE, [var1]) == dna"AGCATGTTAAATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT"

# Test competing mutations (shouldn't happen in reality, but test for weirdness)
@test consensus(REFERENCE, [var1, var2]) ==
    dna"AGCATGTTATATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT"

# Test adjusted threshold
@test consensus(REFERENCE, [var1, var2]; freq=0.9) == REFERENCE

# Test multiple mutations
@test consensus(REFERENCE, [var1, var2, var3]) ==
    dna"AGCATGTTATATAAGATAGGTGTGCTAGTAGGCAGTCAGCGCCAT"
