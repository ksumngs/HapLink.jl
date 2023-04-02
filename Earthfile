VERSION 0.6

docs:
    FROM julia:alpine3.17
    COPY --dir src .
    COPY --dir docs .
    COPY Project.toml .
    RUN apk add --update --no-cache git
    RUN julia --color=yes --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
    RUN julia --color=yes --project=docs -e 'using Documenter: DocMeta, doctest; using HapLink; DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true); doctest(HapLink)'
    RUN julia --color=yes --project=docs docs/make.jl
    SAVE ARTIFACT docs/build AS LOCAL local-output/site

test:
    FOR JULIA_VERSION in 1 1.6 rc
        FROM julia:$JULIA_VERSION
        COPY --dir src .
        COPY --dir docs .
        COPY Project.toml .
        RUN julia --color=yes --project=. -e 'using Pkg; Pkg.instantiate()'
        RUN julia --color=yes --project=. -e 'using Pkg; Pkg.test()'
    END
