# We are going to assume that the raster span the entire globe
longitudes(::Type{T}) where {T <: RasterDataSources.RasterDataSource} = (-180.0, 180.0)
latitudes(::Type{T}) where {T <: RasterDataSources.RasterDataSource} = (-90.0, 90.0)

# CHELSA has a reduced range
longitudes(::Type{CHELSA{BioClim}}) = (-180.0001388888, 179.9998611111)
latitudes(::Type{CHELSA{BioClim}}) = (-90.0001388888, 83.9998611111)

# EarthEnv
longitudes(::Type{T}) where {T <: EarthEnv} = (-180.0, 180.0)
latitudes(::Type{T}) where {T <: EarthEnv} = (-56.0, 90.0)

_left(::Type{T}) where {T <: RasterDataSources.RasterDataSource} = minimum(longitudes(T))
_right(::Type{T}) where {T <: RasterDataSources.RasterDataSource} = maximum(longitudes(T))
_bottom(::Type{T}) where {T <: RasterDataSources.RasterDataSource} = minimum(latitudes(T))
_top(::Type{T}) where {T <: RasterDataSources.RasterDataSource} = maximum(latitudes(T))

# Raster methods
function SimpleSDMPredictor(source::Type{T}, layer; left=_left(T), right=_right(T), bottom=_bottom(T), top=_top(T), kw...) where {T <: RasterDataSources.RasterDataSource}
    rasterpath = RasterDataSources.getraster(source, layer; kw...)
    return _raster(T, rasterpath; left=left, right=right, bottom=bottom, top=top)
end

function _raster(::Type{T}, path::String; bbox...) where {T <: RasterDataSources.RasterDataSource}
    return geotiff(SimpleSDMPredictor, T, path; bbox...)
end

function _raster(::Type{T}, path::Vector{String}; bbox...) where {T <: RasterDataSources.RasterDataSource}
    return [geotiff(SimpleSDMPredictor, T, p; bbox...) for p in path]
end