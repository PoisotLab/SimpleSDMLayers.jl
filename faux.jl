using Revise
using SimpleSDMLayers
using GBIF
using Plots
using Distances
using LinearAlgebra

# One layer
elevation = convert(Float32, SimpleSDMPredictor(WorldClim, Elevation))

# Occurrences
observations = occurrences(
    GBIF.taxon("Hypomyces lactifluorum"; strict=true),
    "hasCoordinate" => "true",
    "country" => "CA",
    "country" => "US",
    "limit" => 300,
)

# Clip raster
elevation = clip(elevation, observations)

# Coordinates
coordinates = [(observations[i].longitude, observations[i].latitude) for i in 1:length(observations)]

_dist = Distances.Haversine()

intras_distances = pairwise(_dist, coordinates)