# `VariationPileup`

```@meta
CurrentModule = HapLink
DocTestSetup = quote
    using HapLink
end
```

```@docs
VariationPileup
```

## Getter methods

```@docs
variation(::VariationPileup)
depth(::VariationPileup)
readpos(::VariationPileup)
quality(::VariationPileup)
strand(::VariationPileup)

```

## Pileup calculations

```@docs
altdepth(::VariationPileup)
frequency(::VariationPileup)
strand_bias(::VariationPileup)
```

## Miscellaneous

```@docs
pileup
```
