# Building an SDM with Bayesian Multivariate Logistic Regression

In this example, we build an SDM using Bayesian Multivariate Logistic Regression (using [Turing.jl]())


```@example mvlogit
using SimpleSDMLayers
using GBIF
using Plots
using LinearAlgebra
using Turing
using StatsFuns: logistic
```

We can get some occurrences for the taxon of interest, _Picea pungens_ or the
Blue Spruce, a conifer native to the Rocky Mountains, although it has been
introduced elsewhere in northeastern North America.

```@example mvlogit
records = occurrences(
    taxon("Carnegiea gigantea"),
    "hasCoordinate" => true,
    "country" => "US"
)

while length(records) < size(records)
    occurrences!(records)
end

scatter(longitudes(records), latitudes(records), lab="")
```

Now we load climate data from worldclim.

```@example mvlogit
function boundingbox(records::GBIFRecords)
    left, right = extrema(longitudes(records)) .+ (-2.0, 2.0)
    bottom, top = extrema(latitudes(records))  .+ (-2.0, 2.0)
    return (left=left, right=right, bottom=bottom, top=top)
end

bounds = boundingbox(records)

predictors = worldclim(collect(1:19); bounds...)

for i in 1:length(environment)
    rescale!(environment[i], (0.,1.))
end

plot(predictors[1])
scatter!(longitudes(records), latitudes(records), lab="")
```

The first thing we want to do is create a `SimpleSDMPredictor` layer that is
$1$ where there is a record of occurrence and $0$ elsewhere.

```@example mvlogit
getOccupancyLayer(envLayer::SimpleSDMPredictor, occupancy) = begin
    latticeLats = latitudes(envLayer)
    latticeLongs = longitudes(envLayer)

    occLayer = similar(envLayer)
    for o in occupancy
        long,lat = longitude(o), latitude(o)
        occLayer[long,lat] = convert(eltype(occLayer), 1)
    end
    occLayer
end

```
Now we need to construct a set of `features` and `labels` to use Flux.

```@example mvlogit
function buildFeaturesMatrix(environment::LT, occurrence::OT) where {LT <: AbstractVector{T} where T <: SimpleSDMPredictor, OT <: GBIFRecords}
    xDim, yDim = size(environment[1])
    numberSpatialPoints = xDim*yDim
    numFeatures = length(environment)
        
    featuresMatrix = zeros(numberSpatialPoints, numFeatures)
    labels = [false for i in 1:numberSpatialPoints]
    occupancyLayer = getOccupancyLayer(environment[1], occurrence)
    
    cursor = 1
    for pt in 1:numberSpatialPoints
        if (!isnothing(environment[1][pt] ))
            for f in 1:numFeatures
                featuresMatrix[cursor, f] = environment[f][pt]  
                labels[cursor] = occupancyLayer[pt]
            end
            cursor += 1
        end
    end
    return (featuresMatrix[1:cursor, :], labels[1:cursor])
end

features, labels = buildFeaturesMatrix(environment, occupancy)
```

Now we define a model in `Flux`, which is as simple as

```@example mvlogit
@model mv_logit(features, labels, σ) = begin
    numDataPts = size(features)[1]
    numFeatures = size(features)[2]
    α ~ Normal(0, σ)
    β ~ MvNormal([0. for i in 1:numFeatures], I*[σ for i in 1:numFeatures])
    for i in 1:numDataPts
        v = α + (β ⋅ features[i,:])
        p = logistic(v)
        labels[i] ~ Bernoulli(p)
    end
end;

```

We can then sample from this model 

```@example mvlogit

chain = sample(mv_logit(features, labels, 1), HMC(0.01, 10), 1500)
plot(chain)
```

and define a function to build a prediction layer based on our its sample

```@example mvlogit

function logit(α, β, features)
    v = α + (β ⋅ features)
    p = logistic(v)
end

function predict(chain, environment)
    predictedProb = similar(environment[1])
    numFeatures = length(environment)

    α = mean(chain[:α])
    βvec = [mean(chain["β[$i]"]) for i in 1:numFeatures]

    for i in 1:length(predictedProb)
        envFeatures= [environment[j][i] for j in 1:numFeatures]

        if (!in(nothing, envFeatures))
            prob::eltype(predictedProb) = (logit(α,βvec, envFeatures))
            predictedProb[i] = prob
        end
    end
    predictedProb
end

predictionLayer = predict(chain, environment)
```


and now we plot

```@example 
pl = plot(predictionLayer, aspectratio=1, size=(700,1000), xlim=(bounds.left, bounds.right))
scatter!(coords, c=:white, alpha=0.2, legend=nothing)

```