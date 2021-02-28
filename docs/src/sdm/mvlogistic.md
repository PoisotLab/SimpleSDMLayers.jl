# Building an SDM with Bayesian Multivariate Logistic Regression

In this example, we build an SDM using Bayesian Multivariate Logistic Regression (using [Turing.jl]())


```@example mvlogit
using SimpleSDMLayers
using GBIF
using StatsPlots
using LinearAlgebra
using Statistics
using Turing
using StatsFuns: logistic
using StatsBase
```

We can get some occurrences for the taxon of interest, _Picea pungens_ or the
Blue Spruce, a conifer native to the Rocky Mountains, although it has been
introduced elsewhere in northeastern North America.

```@example mvlogit
records = occurrences(
    taxon("Picea pungens"),
    "hasCoordinate" => true,
    "country" => "US",
    "decimalLatitude" => (30, 50),
    "limit" => 100
)

while length(records) < min(2000, size(records))
    occurrences!(records)
end

scatter(longitudes(records), latitudes(records), lab="", frame=:box, ratio=1)
```

Now we load climate data from worldclim.

```@example mvlogit
function boundingbox(records::GBIFRecords)
    left, right = extrema(longitudes(records)) .+ (-2.0, 2.0)
    bottom, top = extrema(latitudes(records))  .+ (-2.0, 2.0)
    return (left=left, right=right, bottom=bottom, top=top)
end

predictors = worldclim(1:19; boundingbox(records)...)

for predictor in predictors
    rescale!(predictor, (0.,1.))
end

plot(predictors[1], c=:Greens, frame=:box)
scatter!(longitudes(records), latitudes(records), lab="", msw=0.0, c=:orange, ms=3)
```

The first thing we want to do is create a `SimpleSDMPredictor` layer that is $1$
where there is a record of occurrence and $0$ elsewhere, and 

Now we need to construct a set of `features` and `labels` to use Flux.

```@example mvlogit
labels = Float64.(collect(mask(predictors[1], records)))
features = hcat([collect(predictor) for predictor in predictors]...)

negatives = sample(1:size(features[labels.<1,:],1), 2000, replace=false)
positives = findall(labels.>0)

labels = labels[vcat(negatives, positives)]
features = features[vcat(negatives, positives),:]
```

Now we define a model in `Turing`, which is as simple as

```@example mvlogit
@model mv_logit(features, labels, σ) = begin
    samplesize = size(features,1)
    featurecount = size(features,2)
    α ~ Normal(0.0, σ)
    β ~ MvNormal(fill(0.0, featurecount), I*fill(σ, featurecount))
    for i in 1:samplesize
        v = α + (β ⋅ features[i,:])
        p = logistic(v)
        labels[i] ~ Bernoulli(p)
    end
end;
```

We can then sample from this model 

```@example mvlogit
chain = sample(mv_logit(features, labels, 1.0), HMC(0.01, 10), 1500)
```

and define a function to build a prediction layer based on its sample

```@example mvlogit
logit(α, β, features) = logistic(α + (β ⋅ features))

function predict(chain, predictors)
    prediction = similar(first(predictors))
    featurecount = length(predictors)

    α = mean(chain[:α])
    βvec = [mean(chain["β[$i]"]) for i in 1:featurecount]

    for i in eachindex(prediction)
        local_environment = [predictors[j][i] for j in 1:featurecount]

        if (!in(nothing, local_environment))
            prob::eltype(prediction) = (logit(α, βvec, local_environment))
            prediction[i] = prob
        end
    end
    return prediction
end

prediction = predict(chain, predictors)
```

and now we plot

```@example 
plot(rescale(prediction, collect(0.0:0.05:1.0)), c=:alpine, frame=:box)
scatter!(longitudes(records), latitudes(records), lab="", msw=0.0, c=:black, ms=2, alpha=0.2)
```