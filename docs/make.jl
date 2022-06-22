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
        "API Reference" => [
            "Variation Extensions" => "api/variation.md",
            "VariationInfo" => "api/variationinfo.md",
        ],
    ],
)
deploydocs(; repo="github.com/ksumngs/HapLink.jl", devbranch="master")
