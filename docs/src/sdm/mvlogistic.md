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
Corsican nuthatch, a threatened bird native to Corsica.

```@example mvlogit
records = occurrences(
    taxon("Sitta whiteheadi"),
    "hasCoordinate" => true,
    "limit" => 500
)


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

Now we need to construct a set of `features` and `labels` to use `Turing`.
This creates a `n` by `m` matrix called `features`, where each row corresponds
to a point in the raster, and contains the value of each predictor. We also need
a vector `labels` of length `n`, which correspond to each point in the lattice
and contains `1` if there is an occurrence record and that point, and `0`
otherwise.

```@example mvlogit
labels = Float64.(collect(mask(predictors[1], records)))
features = hcat([collect(predictor) for predictor in predictors]...)

negatives = sample(1:size(features[labels.<1,:],1), 2000, replace=false)
positives = findall(labels.>0)

labels = labels[vcat(negatives, positives)]
features = features[vcat(negatives, positives),:]
```

Now we define `Turing` model to do multivariate logistic regression.

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

Not sure how much in detail I should explain how mvlogit works or assume a base
familiarity and link to other resources.

Effects, represented as the vector $\beta$, are sampled from a multivariate
normal prior with equal variance and no covariance in effects. The predicted
state is computed as

$$\text{logit}(y) = \alpha + \sum_i \Beta_i x_i$$

One could also choose to sample from a set of priors on matrices $\textbf{B}$,
and let the data speak for itself regarding the interaction effects $B_ij$. For
the sake of example we cut out a lot of the fine tuning one would do in an
actual analysis, checking for covariance of predictors and so on.


We can then sample from this model. This line creates a single Markov chain
which runs using Hamilitonian Monte Carlo (`HMC`, see [todo resource on HMC]())
to sample a posterior estimate of our parameters. Here we use `0.01` as the step
size as it resulted in few divergent transitions.

```@example mvlogit
chain = sample(mv_logit(features, labels, 1.0), HMC(0.01, 10), 1500)
```

Now that we have a posterior estimate for the parameters of our model,
we define a function to build a prediction layer based on its sample.

There is almost certaintly a more `SimpleSDMLayers`-y way to do this.

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
            prob::eltype(prediction) = logit(α, βvec, local_environment)
            prediction[i] = prob
        end
    end
    return prediction
end

prediction = predict(chain, predictors)
```

and now we plot our SDM

```@example
plot(rescale(prediction, collect(0.0:0.05:1.0)), c=:alpine, frame=:box)
scatter!(longitudes(records), latitudes(records), lab="", msw=0.0, c=:black, ms=2, alpha=0.2)
```
