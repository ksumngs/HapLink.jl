# Private API Reference

```@meta
CurrentModule = HapLink
DocTestSetup = quote
    using HapLink
end
```

```@autodocs
Modules = [HapLink]
Private = true
Public = false
Filter = f -> f ∉ [
    HapLink.variants,
    HapLink.consensus,
    HapLink.haplotypes,
    HapLink.sequences,
]
```

```@docs
findset
```
