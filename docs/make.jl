using haplink
using Documenter

DocMeta.setdocmeta!(haplink, :DocTestSetup, :(using haplink); recursive=true)

makedocs(;
    modules=[haplink],
    authors="Thomas A. Christensen II <25492070+MillironX@users.noreply.github.com> and contributors",
    repo="https://github.com/MillironX/haplink.jl/blob/{commit}{path}#{line}",
    sitename="haplink.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MillironX.github.io/haplink.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MillironX/haplink.jl",
    devbranch="master",
)
