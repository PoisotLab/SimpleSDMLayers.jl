## Simple layers for Species Distributions Modelling

This package offers very simple types and functions to interact with
bioclimatic data and the output of species distribution models.

[![Build Status](https://travis-ci.org/EcoJulia/SimpleSDMLayers.jl.svg?branch=master)](https://travis-ci.org/EcoJulia/SimpleSDMLayers.jl)

### Installation

For now, installation requires to use the github repo as the package is not tagged yet:

~~~ julia
] add https://github.com/EcoJulia/SimpleSDMLayers.jl#master
~~~

### Type system

All types belong to the abstract `SimpleSDMLayer`, and are organised in the
same way: a `grid` field storing a matrix of data (of any type!), and the
`left`, `right`, `bottom`, and `top` coordinates (as floating point values).

The two core types of the package are `SimpleSDMPredictor` and
`SimpleSDMResponse`. The only difference between the two is that predictors
are immutable, but responses are.

### Methods

Most of the methods are overloads from `Base`. In particular, `SimpleSDMLayer`
objects can be accessed like normal two-dimensional arrays, in which case
they return an object of the same type if called with a range, and the value
if called with a single position.

It is also possible to crop a layer based on a bounding box:

~~~ julia
p[(left, right), (bottom, top)]
~~~

If the layer is of the `SimpleSDMResponse` type, it is possible to write to it:
~~~ julia
p[-74.3, 17.65] = 1.4
~~~

This is only defined for `SimpleSDMResponse`, and `SimpleSDMPredictor`
are immutable.

### Bioclimatic data

#### WorldClim 2.0

The `worldclim` function will get a range, or an array of indices, and return
the corresponding bioclim 2.0 layers at the specified `resolution`. For
example, to get the annual temperature, and annual precipitation:

~~~ julia
temperature, precipitation = worldclim([1,12])
~~~

By default, the function will return the layers for the entire globe, and they
can be cropped later. The layers are returned as `SimpleSDMPredictor` objects.


### Plotting

Using the `Plots` package, one can call the `heatmap`, `density` (requires
`StatsPlots`), and `plot` methods. Note that `plot` defaults to a `heatmap`.

One can also use `scatter(l1, l2)` where both `l1` and `l2` are layers with the
same dimensions and bounding box, to get a scatterplot of the values. Only the
pixels that have non-`NaN` values in *both* layers are shown. Similarly,
`histogram2d` works.
