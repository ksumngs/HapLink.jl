# Contains example data for use in testing
# These strings contain the alignments from section 1.1 "An example" of the SAM v1
# specification, with random ~Q30 quality strings added
module Examples
const SAMStrings = [
    "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t?>@==A?@>AAA??@@?",
    "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t?>@==A?@>AAA??",
    "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t?>@==A?@>AA\tSA:Z:ref,29,-,6H5M,17,0;",
    "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCTTCAGC\t@?>@>@=?=@?",
    "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t?A?@>\tSA:Z:ref,9,+,5S6M,30,1;",
    "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t?A@=>>>@?\tNM:i:1"
]

# These strings are exactly the same as SAMStrings, but they contain an added mutation
# at position 17 in reads 1 and 4
const MutStrings = [
    "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGACACTG\t?>@==A?@>AAA??@@?",
    "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t?>@==A?@>AAA??",
    "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t?>@==A?@>AA\tSA:Z:ref,29,-,6H5M,17,0;",
    "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tACAGCTTCAGC\t@?>@>@=?=@?",
    "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t?A?@>\tSA:Z:ref,9,+,5S6M,30,1;",
    "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t?A@=>>>@?\tNM:i:1"
]
end
