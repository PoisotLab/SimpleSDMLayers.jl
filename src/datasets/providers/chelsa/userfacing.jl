"""
    SimpleSDMPredictor(CHELSA, BioClim, layer::Integer; left=nothing, right=nothing, bottom=nothing, top=nothing)

Download and prepare bioclim layers from the CHELSA database, and returns
them as an array of `SimpleSDMPredictor`s. Layers are called by their number,
from 1 to 19. The list of available layers is given in a table below.
"""
function SimpleSDMPredictor(::Type{CHELSA}, ::Type{BioClim}, layer::Integer=1; kwargs...)
    file = _get_raster(CHELSA, BioClim, layer)
    return geotiff(SimpleSDMPredictor, file; kwargs...)
end

function SimpleSDMPredictor(::Type{CHELSA}, ::Type{BioClim}, mod::CMIP6, fut::SharedSocioeconomicPathway, layer::Integer=1; year="2011-2040", kwargs...)
    file = _get_raster(CHELSA, BioClim, mod, fut, layer, year)
    return geotiff(SimpleSDMPredictor, file; kwargs...)
end
