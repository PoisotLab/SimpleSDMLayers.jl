To facilitate writing julian code, we have overloaded a number of methods from
`Base`. These methods should remove the need to interact with the `grid` field
directly, and also allow to set and get values using the geographic coordinates
(as opposed to the grid positions).

## Layer manipulation

```@docs
convert
copy
collect
eltype
size
stride
similar
Base.show
```

## Access to data and transformations

```@docs
eachindex
getindex
setindex!
replace
replace!
broadcast
```

## Operations on a single layer

```@docs
Base.sum
Base.maximum
Base.minimum
Base.extrema
Base.max
Base.min
Base.sqrt
Base.log
Base.log2
Base.log1p
Base.log10
Base.exp
Base.exp2
Base.exp10
Base.expm1
Base.abs
Statistics.mean
Statistics.median
Statistics.std
Statistics.quantile
```

## Operations on pairs of layers

```@docs
+
-
*
/
Base.hcat
Base.vcat
```
