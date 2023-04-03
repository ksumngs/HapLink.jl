VERSION 0.6
FROM alpine:3.17
COPY --dir src .
COPY --dir docs .
COPY --dir test .
COPY Project.toml .

project:
    COPY --dir src .
    COPY --dir docs .
    COPY --dir test .
    COPY Project.toml .
    SAVE ARTIFACT .

docs:
    FROM julia:alpine3.17
    COPY +project/ .
    RUN apk add --update --no-cache git
    RUN julia --color=yes --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
    RUN julia --color=yes --project=docs -e 'using Documenter: DocMeta, doctest; using HapLink; DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true); doctest(HapLink)'
    RUN julia --color=yes --project=docs docs/make.jl
    SAVE ARTIFACT docs/build AS LOCAL local-output/site

test:
    BUILD +test-latest
    BUILD +test-lts
    BUILD +test-rc

test-latest:
    FROM julia:1
    COPY +project/ .
    RUN julia --color=yes --project=. -e 'using Pkg; Pkg.instantiate()'
    RUN julia --color=yes --project=. -e 'using Pkg; Pkg.test()'

test-lts:
    FROM julia:1.6
    COPY +project/ .
    RUN julia --color=yes --project=. -e 'using Pkg; Pkg.instantiate()'
    RUN julia --color=yes --project=. -e 'using Pkg; Pkg.test()'

test-rc:
    FROM julia:rc
    COPY +project/ .
    RUN julia --color=yes --project=. -e 'using Pkg; Pkg.instantiate()'
    RUN julia --color=yes --project=. -e 'using Pkg; Pkg.test()'
