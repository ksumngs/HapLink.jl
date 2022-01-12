# syntax=docker/dockerfile:1.2

FROM ubuntu:focal

ENV JULIA_VERSION 1.6.5
ENV HAPLINK_VERSION 0.4.3
ENV JULIA_DEPOT_PATH /.julia

# Install the build dependencies
RUN \
  apt-get update && \
  apt-get install -y --no-install-recommends \
    ssh \
    curl \
    git \
    ca-certificates \
    build-essential

# Install Julia
RUN \
  cd / && \
  curl -L "https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" | tar xvz && \
  ln -s /julia-${JULIA_VERSION}/bin/julia /usr/bin/julia

# Install PackageCompiler.jl
RUN \
  julia -e 'using Pkg; Pkg.add("PackageCompiler")'

# Clone and build HapLink.jl
RUN --mount=type=secret,id=SSHKEY \
  cd / && \
  mkdir -p /root/.ssh && \
  ssh-keyscan github.com >> /root/.ssh/known_hosts && \
  cp /run/secrets/SSHKEY /root/.ssh/id_rsa && \
  git clone git@github.com:ksumngs/HapLink.jl.git && \
  rm -rf /root/.ssh && \
  cd HapLink.jl && \
  git checkout v${HAPLINK_VERSION} && \
  julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()' && \
  julia -e 'using PackageCompiler; create_app(".", "build", precompile_execution_file="precompile_app.jl", executables=["haplink" => "haplink"], cpu_target="x86-64")'

FROM ubuntu:focal

COPY --from=0 /HapLink.jl/build/bin /usr/bin
COPY --from=0 /HapLink.jl/build/lib /usr/lib
COPY --from=0 /HapLink.jl/build/share /usr/share

ENTRYPOINT ["/usr/bin/haplink"]
