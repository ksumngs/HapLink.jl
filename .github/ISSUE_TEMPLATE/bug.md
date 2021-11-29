---
name: Bug Report
about: Something isn't working like it's supposed to
labels: bug
---

> _Thanks for helping squash bugs in HapLink!_ :smile:
> _This template is rather extensive. Please fill out all that you can, but
> remember that this template is only a tool: if you feel anything in this
> template is not relevant, simply delete it._

Give a clear and consise description of the bug here

## Expected Behavior

Tell us what you expected to happen

## Current Behavior

Tell us what actually happened

## Steps to Reproduce

Use GitHub's file upload to upload files that can be used to reproduce this bug.

### Input files

- [sample.bam](#)
- [reference.fasta](#)

### Command

```shellsession
haplink sample.bam -r reference.fasta
```

### REPL

```julia
using HapLink
callvariants(myvars...)
```

## Context

Tell us what you are trying to accomplish and/or how this bug affects users in
the "real world."

## Possible Solution

(Optional, but highly useful). If you can figure out a design decision or a
particular function call that is causing the bug, explain here. Even better, if
you know a way to fix it, but don't want to open a pull request, explain that
here.

## Environment

- OS and version: <!-- e.g. CentOS 8, Mac OS (Big Sur), Ubuntu 20.04 -->
- Julia version: <!-- e.g. 1.6.3, 1.6.0, nightly -->
- Julia install source: <!-- e.g. Distro package, Homebrew, Conda, Docker image (include link) -->
- bam-readcount version: <!-- e.g. 1.0.0 -->
- bam-readcount source <!-- e.g. self-built, Conda, Docker image (include link) -->
- HapLink version: <!-- e.g. 0.1.0 -->
- HapLink method: <!-- CLI, REPL, or both -->

Please also locate your `Manifest.toml` file, and link to it here.
