VERSION 0.6
FROM julia:alpine3.17
RUN apk add --update --no-cache git

docs:
    COPY --dir src .
    COPY --dir docs .
    COPY Project.toml .
    RUN julia --color=yes --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
    RUN julia --color=yes --project=docs -e 'using Documenter: DocMeta, doctest; using HapLink; DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true); doctest(HapLink)'
    RUN julia --color=yes --project=docs docs/make.jl
    SAVE ARTIFACT docs/build AS LOCAL local-output/site