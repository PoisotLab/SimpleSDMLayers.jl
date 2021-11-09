# # MultivariateStats.jl integration

# In this example we explore how you can use 
# `MultivariateStats.jl` with `SimpleSDMLayers` 
# to perform [principal-component-analysis (PCA)](wikiling todo) to 
# reduce the dimensionality of a set of `SDMLayers`, or to remove the 
# covariance between layers via data [whitening](linktodo).

# ## Principal-component analysis of many `SDMLayers`

using SimpleSDMLayers
using MultivariateStats
using Plots

# The `SimpleSDMLayers` enables integration with `MultivariateStats.jl` 
# In this example, we will show how this can work.

boundaries =(left = -164.022167, right = -55.250858, bottom = 23.547132, top = 72.105833)
layers = convert.(
    Float32,
    SimpleSDMPredictor(WorldClim, BioClim, 1:19; boundaries...)
)

# fit pca and project to a new set of layers

pca = fit(PCA, layers)
newlayers = transform(pca, layers)

# plot them

pcaplots = plot.(newlayers)
plot(pcaplots...)



# ## Remove correlation between layers (`Whitening`)

# Have to start with two layers that are reasonably correlated.

wlayers = convert.(
    Float32,
    SimpleSDMPredictor(WorldClim, BioClim, 1:2; boundaries...)
)
plot(plot.(wlayers)...)


# Now we call methods just as in `MultivariateStats`

w = fit(Whitening, wlayers)
newlayers = transform(w, wlayers)


# and plot the layers without covar

plot(plot.(newlayers)...)