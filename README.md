## Simple Layers for Species Distributions Modelling

This package offers very simple types and functions to interact with
bioclimatic data and the output of species distribution models.

[![d_stable](https://img.shields.io/badge/Doc-stable-green?style=flat-square)](https://ecojulia.github.io/SimpleSDMLayers.jl/stable/)
[![d_latest](https://img.shields.io/badge/Doc-latest-blue?style=flat-square)](https://ecojulia.github.io/SimpleSDMLayers.jl/latest/)

![version](https://img.shields.io/github/v/tag/EcoJulia/SimpleSDMLayers.jl?sort=semver&style=flat-square)
![CI](https://img.shields.io/github/workflow/status/EcoJulia/SimpleSDMLayers.jl/CI?label=CI&style=flat-square)
![Doc](https://img.shields.io/github/workflow/status/EcoJulia/SimpleSDMLayers.jl/Documentation?label=Doc&style=flat-square)
![Coverage](https://img.shields.io/codecov/c/github/EcoJulia/SimpleSDMLayers.jl?style=flat-square)

[![DOI](https://zenodo.org/badge/187030040.svg)](https://zenodo.org/badge/latestdoi/187030040)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02872/status.svg)](https://doi.org/10.21105/joss.02872)

![gbif_figure](./joss/figures/paper_gbif_1.png)

Curious to know more? Have a look at our [paper in Journal of Open Source Software](https://doi.org/10.21105/joss.02872), our [JuliaCon poster](https://github.com/gabrieldansereau/juliacon-2020-poster/blob/master/juliacon-poster.pdf), our [NextJournal demo notebook](https://nextjournal.com/gabrieldansereau/SimpleSDMLayers-JuliaCon2020-demo/), and our [extended documentation](https://ecojulia.github.io/SimpleSDMLayers.jl/stable/), or keep reading for a quick overview.

### Installation

The currently released version of the package can be installed with:

~~~ julia
] add SimpleSDMLayers
~~~

The package is also designed to work with `GBIF`, so you may want to use the following line instead:

~~~ julia
] add SimpleSDMLayers GBIF
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
p[left=left, right=right, bottom=bottom, top=top]
~~~

If the layer is of the `SimpleSDMResponse` type, it is possible to write to it:
~~~ julia
p[-74.3, 17.65] = 1.4
~~~

This is only defined for `SimpleSDMResponse`, and `SimpleSDMPredictor`
are immutable.

### Bioclimatic data

#### WorldClim 2.1

The `worldclim` function will get a range, or an array of indices, and return
the corresponding bioclim 2.1 layers at the specified `resolution`. For
example, to get the annual temperature, and annual precipitation:

~~~ julia
temperature, precipitation = worldclim([1,12])
~~~

By default, the function will return the layers for the entire globe, and they
can be cropped later. The layers are returned as `SimpleSDMPredictor` objects.


### Plotting

Using the `Plots` package, one can call the `heatmap`, `contour`, `density`
(requires `StatsPlots`), and `plot` methods. Note that `plot` defaults to a
`heatmap`.

~~~ julia
plot(temperature)
~~~

![temperature_figure](./joss/figures/paper_temp_1.png)

One can also use `scatter(l1, l2)` where both `l1` and `l2` are layers with the
same dimensions and bounding box, to get a scatterplot of the values. Only the
pixels that have non-`NaN` values in *both* layers are shown. Similarly,
`histogram2d` works.

## How to contribute

Please read the [Code of Conduct][CoC] and the [contributing guidelines][contr].

[CoC]: https://github.com/EcoJulia/SimpleSDMLayers.jl/blob/master/CODE_OF_CONDUCT.md
[contr]: https://github.com/EcoJulia/SimpleSDMLayers.jl/blob/master/CONTRIBUTING.md
