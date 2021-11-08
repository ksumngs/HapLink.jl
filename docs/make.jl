using haplink
using Documenter

DocMeta.setdocmeta!(haplink, :DocTestSetup, :(using haplink); recursive=true)

makedocs(;
    modules=[haplink],
    authors="Thomas A. Christensen II, Kansas State University, and contributors",
    repo="https://github.com/ksumngs/haplink.jl/blob/{commit}{path}#{line}",
    sitename="haplink.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ksumngs.github.io/haplink.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ksumngs/haplink.jl",
    devbranch="master",
)
