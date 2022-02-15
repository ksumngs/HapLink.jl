# Contains example data for use in testing
# These strings contain the alignments from section 1.1 "An example" of the SAM v1
# specification. Quality strings of 30 were added, with inserted bases having a quality of
# 25 for easy testing
module Examples
const SAMStrings = [
    "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t????????::???????",
    "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t:::??????:????",
    "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t:::::??????\tSA:Z:ref,29,-,6H5M,17,0;",
    "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCTTCAGC\t???????????",
    "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t?????\tSA:Z:ref,9,+,5S6M,30,1;",
    "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t?????????\tNM:i:1",
]

# These strings are exactly the same as SAMStrings, but they contain an added mutation
# at position 17 in reads 1 and 4, and this mutation has a quality score of 35
const MutStrings = [
    "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGACACTG\t????????::??D????",
    "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t:::??????:????",
    "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t:::::??????\tSA:Z:ref,29,-,6H5M,17,0;",
    "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tACAGCTTCAGC\t?D?????????",
    "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t?????\tSA:Z:ref,9,+,5S6M,30,1;",
    "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t?????????\tNM:i:1",
]
end
