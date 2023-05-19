VERSION 0.7
FROM alpine:3.17

docs:
    FROM julia:alpine3.18
    COPY --dir src .
    COPY --dir docs .
    COPY Project.toml .
    RUN apk add --update --no-cache git
    RUN julia --color=yes --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
    RUN julia --color=yes --project=docs -e 'using Documenter: DocMeta, doctest; using HapLink; DocMeta.setdocmeta!(HapLink, :DocTestSetup, :(using HapLink); recursive=true); doctest(HapLink)'
    RUN julia --color=yes --project=docs docs/make.jl
    SAVE ARTIFACT docs/build AS LOCAL local-output/site

test-all:
    FOR JULIA_VERSION IN 1 rc 1.6
        BUILD +test --version $JULIA_VERSION
    END

test:
    ARG version=latest
    ARG project=@.
    ARG precompile=no
    ARG check_bounds=yes
    ARG coverage=true
    ARG depwarn=yes
    ARG force_latest_compatible_version=auto
    ARG inline=yes
    ARG prefix=''
    ARG annotate=false
    ARG runtest_version=1.9.3
    FROM julia:$version
    COPY --dir src .
    COPY --dir docs .
    COPY --dir test .
    COPY Project.toml .
    ENV JULIA_PKG_PRECOMPILE_AUTO=$precompile
    ENV ANNOTATE=$annotate
    ENV COVERAGE=$coverage
    ENV FORCE_LATEST_COMPATIBLE_VERSION=$force_latest_compatible_version
    ENV CHECK_BOUNDS=$check_bounds
    ENV INPUT_DIRECTORIES='src,ext'
    RUN julia --color=yes --project=$project -e ' \
        import Pkg; \
        VERSION >= v"1.5-" && Pkg.Registry.add("General"); \
        VERSION >= v"1.1.0-rc1" ? Pkg.build(verbose=true) : Pkg.build()'
    RUN julia --color=yes -e '\
        if v"1.8pre" < VERSION < v"1.9.0-beta3"; \
            using Pkg; \
            Pkg.activate("tests-logger-env"; shared=true); \
            Pkg.add(Pkg.PackageSpec(name="GitHubActions", version="0.1")); \
        end'
    RUN curl -L https://github.com/julia-actions/julia-runtest/archive/refs/tags/v1.9.3.tar.gz | tar xvz
    RUN julia --color=yes \
        --depwarn=yes \
        --inline=yes \
        --project=@. \
        -e 'include(joinpath("julia-runtest-1.9.3", "test_harness.jl"))'
    RUN curl -L https://github.com/julia-actions/julia-processcoverage/archive/refs/tags/v1.2.2.tar.gz | tar xvz
    RUN julia --color=yes julia-processcoverage-1.2.2/main.jl
    SAVE ARTIFACT lcov.info AS LOCAL local-output/lcov.$version.info
    