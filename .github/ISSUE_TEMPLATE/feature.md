---
name: Feature Request
about: Something isn't working like it's supposed to
labels: enhancement
---

> _Thanks for adding fresh ideas to HapLink!_ :smile:
> _This template is rather extensive. Please fill out all that you can, but
> remember that it is only a tool: if you feel anything in this template is not
> relevant, simply delete it._

Give a clear and consise description of what you want to happen

## Additional Functions

Explain the signature and output of a new API function, in a similar style to
that of Julia docstrings.

## Additional Arguments

Propose a long (`--whole-word`) and short (`-s` single-letter) argument and how
it would affect the analysis and/or output of HapLink.

## Additional Output

Explain which file the output should go in (current output files and
newly-proposed ones are both fair game), what data should be output, and what
format it should be in. If proposing a new format that HapLink doesn't yet
output, please link to the specification for that format.

## Example Usage

Use GitHub's file upload to upload files that can be used as a basis to
construct this feature (if applicable).

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

Tell us what you are trying to accomplish and/or how this feature would affect
users in the "real world."

## Alternatives

Explain what you're doing now, or if any other tools come close to this feature.

## Possible Implementation

(Optional, but highly useful). If you can figure out a good entry point in the
code, or have a one-liner that could be a good starting place, list it here. Even
better, if you know a way to implement it, but don't want to open a pull
request, explain that here.
