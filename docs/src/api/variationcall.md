# `VariationCall`

```@meta
CurrentModule = HapLink
DocTestSetup = quote
    using HapLink
end
```

```@docs
VariationCall
```

## Getter methods

```@docs
variation(::VariationCall)
quality(::VariationCall)
filter(::VariationCall)
depth(::VariationCall)
strand_bias(::VariationCall)
altdepth(::VariationCall)
readpos(::VariationCall)
p_value(::VariationCall)
frequency(::VariationCall)
```

## Variant calling methods

```@docs
variation_test
call_variant
```
