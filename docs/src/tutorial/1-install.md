# [In the beginning](@id install-tutorial)

There are many different ways to install HapLink. Note that some of these
install methods are platform-specific.

```@contents
Pages = ["1-install.md"]
```

## Bioconda

We understand: every bioinformatian is addicted to
[miniconda](https://docs.conda.io/en/latest/miniconda.html). 🐍 And we're not
here to judge. 👩‍⚖️ It's easy and portable and is bundled on most HPCs. If you
already have conda (or [mamba](https://mamba.readthedocs.io/en/latest/)), then
this route is probably for you.

!!! warning
    
    Bioconda install is only supported on Linux

### Install HapLink inside a conda environment

We'll make a new environment with the totally original name "haplink," to house
the new haplink install and add it directly.

```bash
conda create -n haplink -c bioconda -c conda-forge haplink -y
```

### Activate the environment

Now that we've made the environment, let's use it!

```bash
conda activate haplink
```

### Test for errors

Next, cross your fingers 🤞 and run the following command:

```bash
haplink --help
```

The most common error is

```shellsession
The following package could not be installed
└─ haplink   does not exist (perhaps a typo or a missing channel).
```

If this happens,

 1. Check your spelling in the install command
 2. Check that you are using an x86-64 version of conda on Linux

Another common error is

```shellsession
bash: haplink: command not found
bash: /bin/julia: No such file or directory
```

If this happens, check that `CONDA_PREFIX` is set correctly by running
`echo "$CONDA_PREFIX"`, and/or rerun `conda activate haplink`.

If there are no error messages, you're done. You can reuse your
`haplink` environment for the [next tutorial](@ref cli-tutorial).

## Container install

One option for installing HapLink is don't install HapLink. Or rather, pull a
[container](https://apptainer.org/docs/user/1.2/introduction.html#why-use-containers)
that already has HapLink installed, and process files inside of it. HapLink
provides a Docker container that has been tested on [Apptainer](https://apptainer.org).
You should be able to use nearly any container software to run HapLink, but we
recommend Apptainer, due to its ubiquity on HPCs, simple file permissions, and
increased security.

### Download the container

With Apptainer installed, run

```bash
apptainer pull docker://ghcr.io/ksumngs/haplink.jl
```

!!! info "Output"
    
      - haplink.jl_latest.sif

### Run the container as a one-off

You can check to see if the container downloaded correctly by using the
`apptainer exec` command.

```bash
apptainer exec haplink.jl_latest.sif haplink --version
```

### Enter the container to run multiple commands

For more complex commands, it is often better to enter the container's shell
environment and execute commands within the container. Apptainer will include
all files in your working directory as part of the container when doing this.

```shellsession
$ apptainer shell haplink.jl_latest.sif
Apptainer> haplink --version
```

## Julia dependent-install

HapLink is unashamedly a Julia program. If you already have Julia installed,
then you can leverage that existing Julia install to install HapLink.

!!! tip
    
    Under the hood, HapLink can self-install thanks to the power of
    [Comonicon.jl](https://comonicon.org/). Check out their docs if you want to
    learn more, or want to troubleshoot a direct install.

### Check your Julia version

HapLink requires Julia v1.6 or later to run. No exceptions, Kemosabe. Run a
check now to see if you have a high enough version.

```bash
julia --version
```

### Add HapLink to a temporary environment and install

Using the magic of Julia's environments, we can do a "temp install" of the
HapLink package to a temporary directory environment. Because this is a fresh
install, though, it will trigger an installation of the application to a
brand-new isolated environment.

```bash
julia \
  --startup-file=no \
  --history-file=no \
  --quiet \
  -e 'using Pkg; Pkg.activate(;temp=true); Pkg.add(HapLink)'
```

### Add HapLink to PATH

Comonicon will install HapLink, but chances are it won't be accessible yet.
We need to add HapLink to the `PATH`, first.

#### Bash

```bash
echo 'export PATH=$HOME/.julia/bin:$PATH' >> $HOME/.bashrc
. ~/.bashrc
```

#### Zsh

```bash
echo 'export PATH=$HOME/.julia/bin:$PATH' >> $HOME/.zshrc
. ~/.zshrc
```

### Testing

Now, run the `haplink` command.

```bash
haplink --help
```

### Wait, the test failed!

Odds are good that you already had an install of HapLink somewhere in your Julia
depot. If that happens, then you'll have to trigger the install manually. Just
run

```bash
julia \
  --startup-file=no \
  --history-file=no \
  --quiet \
  -e 'using Pkg; Pkg.activate(;temp=true); Pkg.add(HapLink); using HapLink; HapLink.comonicon_install()'
```

then return to the [Add HapLink to PATH](@ref) step.

* * *

Success! We now have a working installation of HapLink. You are now ready to
move on to the [next tutorial](@ref cli-tutorial).
