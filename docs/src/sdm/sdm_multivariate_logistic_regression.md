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

We can get some occurrences for the taxon of interest, _Picea pungens_ or the Blue Spruce,
a conifer native to the Rocky Mountains, although it has been introduced elsewhere in northeastern
North America.


```@example mvlogit
getData(tx::GBIFTaxon; country="US") = begin
    occData = occurrences(tx, "hasCoordinate" => true, "country" => country)
    while (length(occData)< 1000)
        occurrences!(occData)
        @info length(occData)
    end
    occData
end

occurrence = getData(taxon("Picea pungens"))    

plot(occurrence)
```

Now we load climate data from worldclim.

```@example mvlogit

boundingBox(occData) = begin
    left, right = extrema([o.longitude for o in occData]) .+ (-2, 2)
    bottom, top = extrema([o.latitude for o in occData])  .+ (-5, 5)
    (left=left, right=right, bottom=bottom, top=top)
end
bounds = boundingBox(occupancy)

environment = worldclim(collect(1:19); bounds...)
for i in 1:length(environment) rescale!(environment[i], (0.,1.)) end

plot(environment[1])
scatter!(occupancy)
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

chain = mapreduce(c -> sample(mv_logit(features, labels, 1), HMC(0.03, 10), 1000),
    chainscat,
    1
)
end


and define a function to build a prediction layer based on our its sample

```@example mvlogit
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
plot(predictionLayer, frame=:box, aspectratio=1, xlim=(-125,-90))
scatter!(coordinates, c=:white, alpha=0.1, legend=nothing)

```