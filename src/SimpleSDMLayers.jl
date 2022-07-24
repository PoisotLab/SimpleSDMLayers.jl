module SimpleSDMLayers

# This next bit is about being able to change the path for raster assets
# globally, which avoids duplication this argument across multiple functions.
# Important sidenote, this now uses a temporary path, so data will be downloaded
# anew at a new session.
_layers_assets_path = get(ENV, "SDMLAYERS_PATH", tempname())
isdir(_layers_assets_path) || mkpath(_layers_assets_path)

using ArchGDAL
using Downloads
using RecipesBase
using Colors, ColorBlendModes
using ZipFile
using Requires
using Distances
using Statistics
using GeometryBasics
export Point, Polygon
using PolygonOps
using StatsBase
import NetCDF

# Basic types for the package
include(joinpath("lib", "types.jl"))
export SimpleSDMLayer, SimpleSDMResponse, SimpleSDMPredictor

# Main functions to match coordinates
include(joinpath("lib", "coordinateconversion.jl"))

# Implements a series of interfaces (AbstractArray, iteration, and indexing)
include(joinpath("interfaces", "common.jl"))
include(joinpath("interfaces", "iteration.jl"))
include(joinpath("interfaces", "indexing.jl"))
include(joinpath("interfaces", "broadcast.jl"))

# Additional overloads
include(joinpath("lib", "overloads.jl"))

# Raster clipping
include(joinpath("lib", "clip.jl"))

include(joinpath("lib", "generated.jl"))

include(joinpath("lib", "basics.jl"))
export latitudes, longitudes, boundingbox, grid

include(joinpath("datasets", "core.jl"))

include(joinpath("datasets", "readers", "ascii.jl"))
include(joinpath("datasets", "readers", "netcdf.jl"))
include(joinpath("datasets", "readers", "geotiff.jl"))

include(joinpath("datasets", "types.jl"))
export WorldClim, CHELSA, EarthEnv, TerraClimate
export BioClim, LandCover, HabitatHeterogeneity, Elevation, Topography
export PrimaryClimateVariable, SecondaryClimateVariable

# Climate change scenarios
include(joinpath("datasets", "scenarios", "climate.jl"))
export CMIP6, SharedSocioeconomicPathway
export CMIP5, RepresentativeConcentrationPathway
for s in instances(CMIP5)
    @eval export $(Symbol(s))
end
for s in instances(CMIP6)
    @eval export $(Symbol(s))
end
for s in instances(RepresentativeConcentrationPathway)
    @eval export $(Symbol(s))
end
for s in instances(SharedSocioeconomicPathway)
    @eval export $(Symbol(s))
end

include(joinpath("datasets", "layernames.jl"))
export layernames

#include(joinpath("datasets", "chelsa", "download.jl"))
#include(joinpath("datasets", "chelsa", "bioclim.jl"))

#include(joinpath("datasets", "terraclimate", "download.jl"))

# Data interface for WorldClim
include(joinpath("datasets", "providers", "worldclim", "core.jl"))
include(joinpath("datasets", "providers", "worldclim", "download.jl"))
include(joinpath("datasets", "providers", "worldclim", "userfacing.jl"))

# Data interface for EarthEnv
include(joinpath("datasets", "providers", "earthenv", "core.jl"))
include(joinpath("datasets", "providers", "earthenv", "download.jl"))
include(joinpath("datasets", "providers", "earthenv", "userfacing.jl"))

# Data interface for CHELSA
include(joinpath("datasets", "providers", "chelsa", "core.jl"))
include(joinpath("datasets", "providers", "chelsa", "download.jl"))
include(joinpath("datasets", "providers", "chelsa", "userfacing.jl"))

# Pseudoabsences generation
include(joinpath("pseudoabsences", "main.jl"))
include(joinpath("pseudoabsences", "radius.jl"))
include(joinpath("pseudoabsences", "randomselection.jl"))
include(joinpath("pseudoabsences", "surfacerangeenvelope.jl"))
export WithinRadius, RandomSelection, SurfaceRangeEnvelope

include(joinpath("operations", "coarsen.jl"))
include(joinpath("operations", "sliding.jl"))
include(joinpath("operations", "mask.jl"))
include(joinpath("operations", "rescale.jl"))
include(joinpath("operations", "mosaic.jl"))
export coarsen, slidingwindow, mask, rescale!, rescale, mosaic

include(joinpath("recipes", "recipes.jl"))
export bivariate

# Fixes the export of clip when GBIF or others are loaded
export clip

function __init__()
    @require GBIF = "ee291a33-5a6c-5552-a3c8-0f29a1181037" begin
        @info "Loading GBIF support for SimpleSDMLayers.jl"
        include("integrations/GBIF.jl")
    end
    @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" begin
        @info "Loading DataFrames support for SimpleSDMLayers.jl"
        include("integrations/DataFrames.jl")
    end

end

end # module
