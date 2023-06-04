using HapLink
using Documenter

DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true)

makedocs(;
    modules=[HapLink],
    authors="Thomas A. Christensen II, Kansas State University, and contributors",
    repo="https://github.com/ksumngs/HapLink.jl/blob/{commit}{path}#{line}",
    sitename="HapLink.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ksumngs.github.io/HapLink.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "In the beginning (Installation)" => "tutorial/1-install.md",
            "Kicking the tires (Fake sequences)" => "tutorial/2-examples.md",
            "Plays well with others (External tools)" => "tutorial/3-other.md",
            "For advanced beginners (REPL mode)" => "tutorial/4-repl.md",
        ],
        "CLI Reference" => [
            "haplink variants" => "cli/variants.md",
            "haplink consensus" => "cli/consensus.md",
            "haplink haplotypes" => "cli/haplotypes.md",
            "haplink sequences" => "cli/sequences.md",
        ],
        "Public API Reference" => [
            "Variation Extensions" => "api/variation.md",
            "VariationInfo" => "api/variationinfo.md",
            "VariationPileup" => "api/variationpileup.md",
            "VariationCall" => "api/variationcall.md",
            "Consensus Tools" => "api/consensus.md",
            "Pseudoread" => "api/psuedoread.md",
            "Haplotype Extensions" => "api/haplotype.md",
            "HaplotypeCall" => "api/haplotypecall.md",
        ],
        "Private API Reference" => "api/private.md",
    ],
)
deploydocs(; repo="github.com/ksumngs/HapLink.jl", devbranch="master")
