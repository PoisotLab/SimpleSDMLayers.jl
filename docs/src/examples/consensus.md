# Landcover consensus

In this example, we will create a consensus map of landcover for Corsica based
on the EarthEnv data, and measure the variation within each pixel using the
variance. The first step is to load the packages we need, and create a bounding
box:

```@example cons
using SimpleSDMLayers
using Plots

bbox = (left=8.25, right=10.0, bottom=41.2, top=43.2)
```

We will then do two things. First, get the first layer of landcover (see the
help of `landcover` for a list of the layers), and then create a datacube,
organized around dimensions of latitude, longitude, and layer value - we will
only focus on the 11 first variables, since we do not want the information on
open water (layer 12):

```@example cons
bbox = (left=8.25, right=10.0, bottom=41.2, top=43.2)
lc = SimpleSDMPredictor(EarthEnv{LandCover}, 1:12; bbox...)

function safefindmax(x)
    all(isnothing.(x)) && return nothing
    all(isnan.(filter(!isnothing, x))) && return nothing
    replace!(x, nothing => -Inf)
    replace!(x, NaN => -Inf)
    return findmax(x)[2]
end

consensus = convert(Float64, xmosaic(safefindmax, lc, Int64; sanitize=false))

plot(consensus, c=cgrad(:gist_earth, 12, categorical=true))
```

```@example cons
function get_most_common_landuse(f)
    f[isnan.(f)] .= 0.0
    sum(f) == 0 && return NaN
    return last(findmax(f))
end

function shannon(x)
    v = filter(!isnan, x)
    length(v) == 0 && return NaN
    v = v ./ sum(v)
    return -sum(v.*log2.(v))
end
```

```@example cons
bbox = (left=8.25, right=10.0, bottom=41.2, top=43.2)
lc = SimpleSDMPredictor(EarthEnv{LandCover}, 1:12; bbox...)

function safefindmax(x)
    all(isnothing.(x)) && return nothing
    all(isnan.(filter(!isnothing, x))) && return nothing
    replace!(x, nothing => -Inf)
    replace!(x, NaN => -Inf)
    return findmax(x)[2]
end

function safeshannon(x)
    all(isnothing.(x)) && return nothing
    all(isnan.(filter(!isnothing, x))) && return nothing
    replace!(x, nothing => 0.0)
    replace!(x, NaN => 0.0)
    filter!(!isequal(0.0), x)
    sum(x) == 0.0 && return nothing
    length(x) == 0 && return nothing
    x = x ./ sum(x)
    return -sum(x.*log2.(x))
end

consensus = convert(Float64, xmosaic(safefindmax, lc, Int64; sanitize=false))
het = convert(Float64, xmosaic(safeshannon, lc, Float64; sanitize=false))
replace!(het, -0.0 => nothing)

plot(het)

plot(consensus, c=cgrad(:gist_earth, 12, categorical=true))
```

```@example cons
p1 = plot(consensus, c=cgrad(:Set3_11, categorical=true), frame=:none)
p2 = plot(entropy, c=:bamako, frame=:none)

plot(p1, p2, size=(900, 400), dpi=600)
```
