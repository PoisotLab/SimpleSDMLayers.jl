module SimpleSDMLayers

using GDAL
using HTTP
using RecipesBase
using ZipFile
using Requires

include(joinpath("lib", "types.jl"))
export SimpleSDMLayer, SimpleSDMResponse, SimpleSDMPredictor

include(joinpath("lib", "overloads.jl"))

include(joinpath("lib", "generated.jl"))

include(joinpath("lib", "basics.jl"))
export latitudes, longitudes

include(joinpath("lib", "geotiff.jl"))
export geotiff

include(joinpath("bioclimaticdata", "worldclim.jl"))
export worldclim

include(joinpath("operations", "coarsen.jl"))
export coarsen

include(joinpath("recipes", "recipes.jl"))

function __init__()
    @require GBIF="ee291a33-5a6c-5552-a3c8-0f29a1181037" begin
        @info "GBIF integration loaded"
        include(joinpath(dirname(pathof(SimpleSDMLayers)), "integrations", "GBIF.jl"))
    end
end

end # module
