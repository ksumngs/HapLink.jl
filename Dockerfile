# syntax=docker/dockerfile:1.2

FROM mgibio/bam-readcount:1.0.0

ENV JULIA_VERSION 1.6.4

RUN --mount=type=secret,id=SSHKEY \
  apt-get update && \
  apt-get install -y --no-install-recommends \
    ssh \
    curl \
    git \
    ca-certificates \
    build-essential && \
  mkdir -p /tmp/julia && \
  curl -L "https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" | tar xvz -C /tmp/julia && \
  /tmp/julia/julia-${JULIA_VERSION}/bin/julia -e 'using Pkg; Pkg.add("PackageCompiler")' && \
  mkdir -p /root/.ssh && \
  ssh-keyscan github.com >> /root/.ssh/known_hosts && \
  cp /run/secrets/SSHKEY /root/.ssh/id_rsa && \
  git clone git@github.com:ksumngs/HapLink.jl.git && \
  rm -rf /root/.ssh && \
  cd HapLink.jl && \
  git checkout v0.2.0 && \
  /tmp/julia/julia-${JULIA_VERSION}/bin/julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()' && \
  /tmp/julia/julia-${JULIA_VERSION}/bin/julia -e 'using PackageCompiler; create_app(".", "build", precompile_execution_file="precompile_app.jl", executables=["haplink" => "haplink", "make-haplotype-fastas" => "make_haplotype_fastas"], cpu_target="x86-64")' && \
  cp -r build/bin/* /usr/bin && \
  cp -r build/lib/* /usr/lib && \
  cp -r build/share/* /usr/share && \
  cd .. && \
  rm -rf /tmp/julia HapLink.jl

ENTRYPOINT ["/usr/bin/haplink"]
