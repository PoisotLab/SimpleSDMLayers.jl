"""
    LayerProvider

A `LayerProvider` is an abstract type used to dispatch the correct call of
`SimpleSDMPredictor` to a specific dataset. A dataset is specified by a
`LayerProvider` and a `LayerDataset`, as well as optionally one or multiple
layers, and future climate information, resolution, or dates.
"""
abstract type LayerProvider end

"""
    LayerDataset

A `LayerDataset` is a specific set of rasters provided by a `LayerProvider`. For
a number of dataset types that are very broad (`LandCover`,
`HabitatHeterogeneity`), the precise mapping of layers is documented in their
`SimpleSDMPredictor` method.
"""
abstract type LayerDataset end

_rasterpath(::Type{<:LayerProvider}) = "Providers"
_rasterpath(::Type{<:LayerDataset}) = "Datasets"
_rasterpath(pr::Type{<:LayerProvider}, ds::Type{<:LayerDataset}) = joinpath(_rasterpath(pr), _rasterpath(ds))

"""
    provides(::Type{<:LayerProvider}, ::Type{<:LayerDataset})

This function *must* be overloaded for each provider/dataset pair, to return
true for the valid combinations.
"""
function provides(::Type{<:LayerProvider}, ::Type{<:LayerDataset})
    return false
end

"""
    provides()

Returns the list of datasets provided by each know data provider.
"""
function provides()

end

function layernames(::Type{<:LayerProvider}, ::Type{<:LayerDataset})
    return ()
end

layernames(p::Type{<:LayerProvider}, d::Type{<:LayerDataset}, i::Int) = layernames(p, d)[i]
layernames(p::Type{<:LayerProvider}, d::Type{<:LayerDataset}, i::AbstractArray) = layernames(p, d)[i]

# Default function to read a layer is geotiff
_readfunction(::Type{<:LayerProvider}, ::Type{<:LayerDataset}) = SimpleSDMLayers.geotiff
