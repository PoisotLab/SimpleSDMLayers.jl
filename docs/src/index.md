# Simple SDM Layers in *Julia*

The `SimpleSDMLayers` provides an interface to facilitate the manipulation of
raster data for species distributions modeling (and possibly other applications)
in *Julia*.

The two core types of the package are `SimpleSDMPredictor` and
`SimpleSDMResponse`. The only difference between the two is that predictors are
immutable, but responses are. All types belong to the abstract `SimpleSDMLayer`,
and are organised in the same way: a `grid` field storing a matrix of data (of
any type!), and the `left`, `right`, `bottom`, and `top` coordinates (as
floating point values).

Note that both types are parametric, *i.e.* `SimpleSDMPredictor{Float32}` has
`Float32`-valued cells. Internally, cells that do not have a value are
represented as `nothing`. This allows the type of the layer to be whatever the
user needs - numbers, symbols, or other types.