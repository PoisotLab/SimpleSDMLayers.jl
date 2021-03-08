# Landcover consensus

In this example, we will create a consensus map of landcover for Corsica based
on the EarthEnv data, and measure the variation within each pixel using
Shannon's entropy. We will use the `mosaic` function, which is designed
specifically for these things.

```@example cons
ENV["RASTERDATASOURCES_PATH"] = "rasterdata"
using SimpleSDMLayers
using Plots
default(; frame=:box)

bbox = (left=8.25, right=10.0, bottom=41.2, top=43.2)
lc = SimpleSDMPredictor(EarthEnv{LandCover}, 1:12; bbox...)
```

The first step is to measure Shannon's entropy, which we can get by only
counting the values that are not `nothing`:

```@example cons
function shannon(x)
    x = x ./ sum(x)
    return -sum(x.*log2.(x))
end

lc_heterogeneity = mosaic(shannon, lc, Float32; sanitize=true)
```

The `sanitize=true` keyword (`true` is the default) is here to make sure that
all of the values that are `NaN`, `nothing`, etc, in the input layers are *not*
passed to the function.

We can visualize the result of this operation:

```@example cons
plot(lc_heterogeneity, c=:dense, xlab="Longitude", yab="Latitude")
```

In the next step, we want to return the identifier of the layer with the highest
value, to get the dominant landcover type. For this operation, we need *all*
values includin `nothing` and `NaN`, so we will write a "safe" `findmax`
function:

```@example cons
function safefindmax(x)
    all(isnothing.(x)) && return nothing
    all(isnan.(filter(!isnothing, x))) && return nothing
    replace!(x, nothing => -Inf)
    replace!(x, NaN => -Inf)
    return findmax(x)[2]
end

lc_palette = [colorant"#32a852", colorant"#4cc76c", colorant"#4c702a", colorant"#1c8054", colorant"#a2ad58", colorant"#c2b057", colorant"#6e5f16", colorant"#227375", colorant"#751010", colorant"#c1f6f7", colorant"#999999", colorant"#1226ff"]

# We get the values as Float64 for plotting
lc_consensus = mosaic(safefindmax, lc, Float64; sanitize=false)

plot(lc_consensus, c=cgrad(lc_palette, [1, 12], categorical=true))
xaxis!("Longitude")
yaxis!("Latitude")
```