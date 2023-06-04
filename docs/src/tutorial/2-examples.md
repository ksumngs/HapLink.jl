# [Kicking the tires](@id cli-tutorial)

At this point, we'll play with the example sequences included _gratis_ 💰 with
HapLink. No, they don't represent anything ☁️, and they aren't particularly
interesting 🥱, but they **do** run fast 🏇, so we can get a handle on how the
interface and workflow operate.

```@contents
Pages = ["2-examples.md"]
```

## Getting the goods

Let's get the example files from the code repository. In your terminal, run

```bash
wget https://github.com/ksumngs/HapLink.jl/raw/v1.0.0-rc1/example/reference.fasta
wget https://github.com/ksumngs/HapLink.jl/raw/v1.0.0-rc1/example/sample.bam
wget https://github.com/ksumngs/HapLink.jl/raw/v1.0.0-rc1/example/sample.bam.bai
```

!!! info "Output"
    
      - reference.fasta
      - sample.bam
      - sample.bam.bai

## Spot the difference

In order for HapLink to call haplotypes, it needs to know which sequence
differences are due to sequencing errors, and which are due to genetic mutation.
This process is known as _variant calling_, and HapLink comes bundled with a
variant caller for just this type of occasion, which requires the reference
genome and the alignment BAM file. Since we have both of those, let's run
variant calling now.

```bash
haplink variants reference.fasta sample.bam
```

!!! info "Output"
    
    _None_

HapLink by default outputs to standard output, so the variant calls were printed
on your screen instead of saved 😡. That's okay, though 😌. It's often good to
visually check your variant calls, and it this case we absolutely needed to.
Notice that none of the variants got a `PASS` filter. In fact, all of them were
weeded out by too high of thresholds for depth (remember we only have 10
sequences) and significance. Let's readjust (and save our results this time).

```bash
haplink \
  variants \
  reference.fasta \
  sample.bam \
  --depth 3 \
  --significance 0.1 \
  | tee sample.vcf
```

!!! info "Output"
    
      - sample.vcf

These settings seemed to work out well. Let's stick with them and move on.

## The general lay of the land

At this point, we're going to take a break from haplotype calling and convert
those variant calls into a useful summary: the consensus sequence. HapLink can
call the consensus sequence based solely off variant calls. Let's see that now.

```bash
haplink consensus reference.fasta sample.vcf | tee sample.consensus.fasta
```

!!! info "Output"
    
      - sample.consensus.fasta

## The star attraction

And now it's time for haplotype calling. Before you get your hopes up, there are
no _true_ haplotypes in this file. If 10 reads could yield subconsenus
mysteries, then bioinformatics would be a super easy job. Alas, we live in the
real world, and we'll have to stretch mathematical constructs to get anything
out of these reads.

```bash
haplink \
  haplotypes \
  reference.fasta \
  sample.vcf \
  sample.bam \
  --consensus-frequency 0.75 \
  | tee sample.yml
```

!!! info "Output"
    
      - sample.yml

You can see that HapLink found only one haplotype in this alignment, but
(spoilers!) this isn't really a haplotype. This is just the consensus sequence,
formatted in HapLink's haplotype scheme. The first haplotype in any output file
is always the consensus sequence.

## Haplotypes in the Matrix

If you have reads that don't span the entire genome (like we have here), you can
use HapLink's maximum likelihood simulator to "create" full-length reads by
splicing together reads and look for haplotypes on them. Even though we know
there aren't any haplotypes in this sample, let's get out the simulator and give
it a try.

```bash
haplink \
  haplotypes \
  reference.fasta \
  sample.vcf \
  sample.bam \
  --consensus-frequency 0.75 \
  --simulated-reads \
  --iterations 100 \
  --overlap-min 0 \
  --overlap-max 100 \
  | tee sample.ml.yml
```

!!! info "Output"
    
      - sample.ml.yml

Still nothing, huh? Like I said, no haplotypes here, and simulation can't change
that. Note that simulating full-length reads used _a lot_ more computational
power, so you should try to stick with full-length reads when you can!

## But, what does it mean?

HapLink's haplotype YAML files contain everything needed to recreate the
haplotype computation, but they can't really be used by any other programs.
That's why there's the `sequences` command, so haplotype sequences can be saved
into FASTA format for use by other tools. Let's try this now.

```bash
haplink sequences reference.fasta sample.yml > sample.fasta
```

!!! info "Output"
    
      - sample.fasta

The output only contains the consensus sequence. Since we knew that there were
no haplotypes in this sample we could have used `haplink consensus`, instead to
get the same result.

* * *

You are now ready to move on to the [next tutorial](@ref integration-tutorial)!
