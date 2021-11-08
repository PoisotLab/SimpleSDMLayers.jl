# # Principal-component analysis of many `SDMLayers`

using SimpleSDMLayers
using MultivariateStats
using Plots


# The `SimpleSDMLayers` enables integration with `MultivariateStats.jl` 
# In this example, we will show how this can work.

boundaries = (left=-12.0, right=30.0, bottom=36.0, top=72.0)

layers = convert.(
    Float32,
    SimpleSDMPredictor(WorldClim, BioClim, 1:19; boundaries...),
)

plot(layers[10])

pca = fit(PCA, layers)
newlayers = transform(pca, layers)



plot(plot.(newlayers, size=(300,300))...)
