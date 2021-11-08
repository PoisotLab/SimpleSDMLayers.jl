# # Principal-component analysis of many `SDMLayers`

using SimpleSDMLayers
using MultivariateStats
using Plots

# The `SimpleSDMLayers` enables integration with `MultivariateStats.jl` 
# In this example, we will show how this can work.

boundaries =(left = -164.022167, right = -55.250858, bottom = 23.547132, top = 72.105833)

layers = convert.(
    Float32,
    SimpleSDMPredictor(WorldClim, BioClim, 1:19; boundaries...),
)


pca = fit(PCA, layers)
newlayers = transform(pca, layers)

# idk why this doesn't work

pcaplots = plot.(newlayers)
plot(pcaplots...)
