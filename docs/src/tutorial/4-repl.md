# [For advanced beginners (REPL mode)](@id repl-tutorial)

Julia is an ahead-of-time compiled language. Practically, that means that every
time you restart Julia, you have to recompile all the code you were running.
Using HapLink on the command line involves up to four different commands.
Translation: up to four cases where you lose time to recompiling code that was
just running. Surely there's a better way, right? Well, you can stay within a
single Julia session by using HapLink's REPL mode.

!!! tip
    
    Julia's latency (aka, Time-to-first-plot or TTFP) is a big deal among Julia
    programmers. Although there's no "definitive" place to learn about TTFP,
    [Jakob Nissen's blog](https://viralinstruction.com/posts/badjulia/#compile_time_latency)
    provides some great explanations and actionable advice for reducing latency.

```@contents
Pages = ["4-repl.md"]
```

## Installing HapLink ... into a project

Wait, isn't that what [the first tutorial](@ref install-tutorial) was about?
Yes and no. Those projects installed HapLink's command line interface to a
global environment. We're going to locally install HapLink as a library into a
local package environment.

!!! tip
    
    If you want to learn more about package environments, check out the
    [official Julia environments documentation](https://pkgdocs.julialang.org/v1/environments/)
    or [Julius Krumbiegel's tutorial](https://jkrumbiegel.com/pages/2022-08-26-pkg-introduction/index.html).

Press `]` to enter Pkg mode, then follow along with these prompts.

```plaintext
(@v1.6) pkg> activate .
  Activating project at `~/haplink-tutorial-4`

(haplink-tutorial-4) pkg> add HapLink
```

## Finding the example files

Instead of downloading the example files from the internet, we can reference
them directly from the package.

```@repl main
using HapLink
const PACKAGE_DIR = dirname(dirname(pathof(HapLink)))
const EXAMPLE_DIR = joinpath(PACKAGE_DIR, "example")
const EXAMPLE_BAM = joinpath(EXAMPLE_DIR, "sample.bam")
const EXAMPLE_REFERENCE = joinpath(EXAMPLE_DIR, "reference.fasta")
```

Perfect!

## Pileup time

Variants can't exist as a single unit, they have to be supported by a set of
reads. We can view those reads in a _pileup_, which can easily be made by using
the [`pileup`](@ref) function.

```@repl main
sample_pileup = pileup(EXAMPLE_BAM, EXAMPLE_REFERENCE)
```

## But is it real?

Pileups tell us what the sequences are saying, but they can't tell us if there
are any real variants. We can use the [`call_variant`](@ref) function to filter
a pileup into a set of variants.

```@repl main
all_variant_calls = VariationCall[]
for pile in sample_pileup
    push!(all_variant_calls, call_variant(pile, 0.5; D=0x02, Q=10.0, X=0.01, F=0.05))
end
all_variant_calls
```

We can now remove all variants that did not pass our filters.

```@repl main
filter!(var -> filters(var) == ["PASS"], all_variant_calls)
```

## Rule of the majority

### Filtering variants

We can sort into consensus and subconsensus variants. Note that for this
fake set, consensus variants will be taken at a >75% frequency.

```@repl main
consensus_variant_calls = filter(var -> frequency(var) > 0.75, all_variant_calls)
subconsensus_variant_calls = filter(
    var -> !(var in consensus_variant_calls), all_variant_calls
)
consensus_variations = variation.(consensus_variant_calls)
subconsensus_variations = variation.(subconsensus_variant_calls)
```

### Creating the consensus sequence

We will now need to convert the consensus variations into a consensus sequence.
HapLink doesn't have these tools in-and-of itself, so we'll have to reach into
its dependency packages.

```plaintext
(haplink-tutorial-4) pkg> add BioAlignments SequenceVariation
   Resolving package versions...
    Updating `~/haplink-tutorial-4/Project.toml`
  [00701ae9] + BioAlignments v3.1.0
  [eef6e190] + SequenceVariation v0.2.2
  No Changes to `~/haplink-tutorial-4/Manifest.toml`
```

```@repl main
using BioAlignments, SequenceVariation
reference_sequence = reference(first(consensus_variations))
consensus_haplotype = Haplotype(reference_sequence, consensus_variations)
consensus_sequence = reconstruct(consensus_haplotype)
consensus_alignment = PairwiseAlignment(
    AlignedSequence(consensus_sequence, Alignment(HapLink.cigar(consensus_haplotype))),
    reference_sequence,
)
```

### Minority report

We are interested in looking at the subconsensus variants, but they need to be
reoriented to be put in terms of the consensus sequence.

```@repl main
map!(
    var -> translate(var, consensus_alignment),
    subconsensus_variations,
    subconsensus_variations,
)
```

## Bringing in the reads

Now that we have a consensus sequence, we can properly import the reads for
haplotype calling into HapLink's specialized [`Pseudoread`](@ref) class.
There is a convenient [`pseudoreads`](@ref) function that can directly convert a
BAM file for us.

```@repl main
haplotyping_pseudoreads = pseudoreads(EXAMPLE_BAM, consensus_sequence)
```

### Cutting the chaff

We only need reads that span every one of our subconsensus variants, so let's
get rid of every other read.

```plaintext
(haplink-tutorial-4) pkg> add BioGenerics
   Resolving package versions...
    Updating `~/haplink-tutorial-4/Project.toml`
  [47718e42] + BioGenerics v0.1.2
  No Changes to `~/haplink-tutorial-4/Manifest.toml`
```

```@repl main
using BioGenerics
filter!(
    var -> leftposition(var) >= leftposition(first(subconsensus_variations)),
    haplotyping_pseudoreads,
)
filter!(
    var -> rightposition(var) >= rightposition(first(subconsensus_variations)),
    haplotyping_pseudoreads,
)
```

## Finding haplotypes...

This one is a bit of a mind-bender, but bear with me.

```@repl main
haplotyping_haploypes = haplotype.(haplotyping_pseudoreads)
valid_haplotypes = findset(
    subconsensus_variations, hap -> ishaplotype(hap, haplotyping_haploypes)
)
```

## ...and deciding if they're real

Just because a haplotype _could_ exist, doesn't mean it should. You can check
the statistics of a haplotype by converting it into a [`HaplotypeCall`](@ref).

```@repl main
map(hap -> HaplotypeCall(hap, haplotyping_haploypes), valid_haplotypes);
hc = first(ans)
depth(hc)
linkage(hc)
significance(hc)
```

Nope! This one sure isn't!

* * *

This tutorial only scratches the surface of what's possible with HapLink's REPL
mode. Take a look at the full API reference for more details on how to tweak
your haplotype calling. Have fun!
