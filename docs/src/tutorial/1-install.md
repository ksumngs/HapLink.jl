# [In the beginning](@id install-tutorial)

There are many different ways to install HapLink. Here we walk you through two
of the most common. If you're one of the 0.01% who needs a different method,
then we trust you can extrapolate from these instructions. Note that all of
these tutorials assume you have a Unix-type system (MacOS, BSD, Linux). Windows
command-line support is basically non-existant!

```@contents
Pages = ["1-install.md"]
```

## Bioconda

We understand: every bioinformatian is addicted to
[miniconda](https://docs.conda.io/en/latest/miniconda.html). 🐍 And we're not
here to judge. 👩‍⚖️ It's easy and portable and is bundled on most HPCs. If you
already have conda (or [mamba](https://mamba.readthedocs.io/en/latest/)), then
this route is probably for you.

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

Check for error messages, but otherwise you're done. You can reuse your
`haplink` environment for the [next tutorial](@ref cli-tutorial).

## Comonicon

HapLink is unashamedly a Julia program. If you already have Julia installed,
then you can leverage that existing Julia install to install HapLink thanks to
the power of [Comonicon.jl](https://comonicon.org/).

### Check your Julia version

HapLink requires Julia v1.6 or later to run. No exceptions, Kemosabe. Run a
check now to see if you have a high enough version.

```bash
julia --version
```

### Add HapLink to a temporary environment and install

Using the magic 🪄 of Julia's environments, we can do a "temp install" of the
HapLink package to a temporary directory environment. Because this is a fresh
install, though, it will trigger Comonicon to install the application to a
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
