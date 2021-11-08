# # Principal-component analysis of many `SDMLayers`

# The `SimpleSDMLayers` enables integration with `MultivariateStats.jl` 
# In this example, we will show how this can work.

boundaries = (left=-12.0, right=30.0, bottom=36.0, top=72.0)

layers = convert(
    Float16,
    SimpleSDMPredictor(WorldClim, HabitatHeterogeneity, 1:19; resolution=5, boundaries...),
)

using MultivariateStats

pca = fit(PCA, layers)
newlayers = transform(pca, layers)