# Writing BIOCLIM from scratch

In this example, we will write the BIOCLIM species distribution model using
`SimpleSDMLayers.jl` and `GBIF.jl`.

```@example bioclim
using SimpleSDMLayers
using GBIF
using Plots
using StatsBase
using Statistics
```

We can get some occurrences for the taxon of interest:

```@example bioclim
obs = occurrences(
    GBIF.taxon("Alces alces", strict=true),
    "hasCoordinate" => "true",
    "continent" => "EUROPE",
    "limit" => 50
)
while length(obs) < 500
    occurrences!(obs)
end
```

This query uses a range for the longitude and latitude, so as to make sure that
we get a relatively small region. Before we get the layers, we will figure out
the bounding box for the observations - just to make sure that we will have
something large enough, we will add a 2 degrees padding around it:

```@example bioclim
left, right = extrema([o.longitude for o in obs]) .+ (-5,5)
bottom, top = extrema([o.latitude for o in obs]) .+ (-5,5)
```

With this information in hand, we can start getting our variables. In this
example, we will take all worldclim data, at the default 10 arc minute
resolution:

```@example bioclim
predictors = worldclim(1:19; left=left, right=right, bottom=bottom, top=top);
first(predictors)
```

The point of BIOCLIM (the model, not the dataset) is that the score assigned to
a pixel is maximal is this pixel is the *median* value for a given variable -
therefore, we need to measure the cumulative density function for every pixel in
every variable.

```@example bioclim
_pixel_score(x) = 2.0(x > 0.5 ? 1.0-x : x)

function SDM(layer::T, observations::GBIFRecords) where {T <: SimpleSDMLayer}
    qf = ecdf(layer[observations]) # We only want the observed values
    return (_pixel_score∘qf)
end
```

Note that we use the ∘ (`\circ`) operator to chain the quantile estimation and
the pixel scoring, which requires Julia 1.4. This function returns a *model*,
*i.e.* a function that we can broadcast to a given layer, which might not be the
same one we used for the training.

The next step in BIOCLIM is to get the *minimum* suitability across all layers
for every pixel. Because we have a `min` method defined for a pair of layers, we
can call `minimum` on an array of layers:

```@example bioclim
function SDM(predictors::Vector{T}, models) where {T <: SimpleSDMLayer}
    @assert length(models) == length(predictors)
    return minimum([broadcast(models[i], predictors[i]) for i in 1:length(predictors)])
end
```

The advantage of this approach is that we can call the `SDM` method for
prediction on a smaller layer, or a different layer. This can allow us to do
thing like stitching layers together with `hcat` and `vcat` to use
multi-threading, or use a different resolution for the prediction than we did
for the training.

```@example bioclim
models = [SDM(predictor, obs) for predictor in predictors]
prediction = SDM(predictors, models)
```

Just because we may want to visualize this result in a transformed way, *i.e.*
by looking at the quantiles of suitability, we can call the `rescale!` function:

```@example bioclim
rescale!(prediction, collect(0.0:0.01:1.0))
```

As this map now represents the quantiles of suitability, we may want to remove
the lower 5%. For this, we need to create a boolean mask, which we can do by
broadcasting a conditional:

```@example bioclim
cutoff = broadcast(x -> x > 0.05, prediction)
```

This map can be plotted as we would normally do:

```@example bioclim
plot(prediction, frame=:box, c=:lightgrey) # Plot a uniform background
plot!(mask(cutoff, prediction), clim=(0,1), c=:bamako)
scatter!([(o.longitude, o.latitude) for o in obs], ms=4, c=:orange, lab="")
xaxis!("Longitude")
yaxis!("Latitude")
```

And there it is! A simple way to write the BIOCLIM model by building on the
integration between SimpleSDMLayers and GBIF.