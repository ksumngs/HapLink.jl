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

@test consensus_variants([var1, var2, var3]) == [var2, var3]
@test consensus_variants([var1, var2, var3]; freq=0.69) == [var2]
