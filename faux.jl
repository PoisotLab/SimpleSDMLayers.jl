using Revise
using SimpleSDMLayers
using GBIF
using Plots
using Distances
using LinearAlgebra
using StatsBase

# One layer - Prince Edward Island, so the geography is actually complex
_bbox = (left=-64.5, right=-61.9, bottom=45.9, top=47.1)
elevation = convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; _bbox..., resolution=0.5))
plot(elevation, frame=:box, c=:bamako, dpi=400)

# Occurrences
observations = occurrences(
    GBIF.taxon("Vulpes vulpes"; strict=true),
    "hasCoordinate" => "true",
    "decimalLatitude" => (45.9, 47.1),
    "decimalLongitude" => (-64.5, -61.9),
    "limit" => 300,
)

"""
Get the coordinates for a list of observations, filtering the ones that do not
correspond to valid layer positions
"""
function coordinates(observations, layer)
    xy = [(observations[i].longitude, observations[i].latitude) for i in 1:length(observations)]
    filter!(c -> !isnothing(layer[c...]), xy)
    return xy
end

"""
Get the distance between points as a matrix
"""
function distance_matrix(xy; D=Distances.Haversine(6371.0))
    return pairwise(D, xy)
end

"""
Bin a distance matrix, where m is the maximum distance allowed
"""
function bin_distances(D, m)
    b = fit(Histogram, vec(D) ./ m, 0.0:0.05:1.0).weights
    return b ./ sum(b)
end

"""
This solves the direct (first) geodetic problem assuming Haversine distances are
a correct approximation of the distance between points.
"""
function randompoint(ref, d; R=6371.0)
    # Convert the coordinates from degrees to radians
    λ, φ = deg2rad.(ref)
    # Get the angular distance
    δ = d / R
    # Pick a random bearing (angle w.r.t. true North)
    α = deg2rad(rand() * 360.0)
    # Get the new latitude
    φ2 = asin(sin(φ) * cos(δ) + cos(φ) * sin(δ) * cos(α))
    # Get the new longitude
    λ2 = λ + atan(sin(α) * sin(δ) * cos(φ), cos(δ) - sin(φ) * sin(φ2))
    # Return the coordinates in degree
    return rad2deg.((λ2, φ2))
end

"""
Get a random point given a layer and the observed distance matrix
"""
function initial_point(layer, D; R=6371.0)
    invalid = true
    global random_destination
    while invalid
        # Get a random distance
        random_distance = rand(D)
        random_starting_point = rand(keys(layer))
        random_destination = randompoint(random_starting_point, random_distance; R=R)
        invalid = isnothing(layer[random_destination...])
    end
    return random_destination
end

layer = copy(elevation)
xy = coordinates(observations, layer)
D = distance_matrix(xy)
m = maximum(D)
b = bin_distances(D, m)

"""
Generates the initial proposition for points
"""
function _initial_proposition(layer, xy, D)
    i0 = initial_point(layer, D)
    all_points = fill(i0, length(xy))
    for i in 2:length(xy)
        global proposition
        invalid = true
        while invalid
            proposition = initial_point(layer, D)
            invalid = isnothing(layer[proposition...]) || (maximum(distance_matrix(all_points[1:i])) > maximum(D))
        end
        all_points[i] = proposition
    end
    return all_points
end

mocks = _initial_proposition(layer, xy, D)
binned_distances = bin_distances(D, maximum(D))

function _points_distance(mocks, m, b)
    mD = distance_matrix(mocks)
    maximum(mD) > m && return Inf
    mb = bin_distances(mD, m)
    return kl_divergence(mb, b)
end

function _improve_one_point!(mocks, layer, D, binned_distances, d0)
    random_distance = rand(D)
    random_point = rand(eachindex(mocks))
    global proposition, nd
    invalid = true
    while invalid
        proposition = randompoint(mocks[random_point], random_distance)
        invalid = isnothing(layer[proposition...])
    end
    current = mocks[random_point]
    mocks[random_point] = proposition
    dt = _points_distance(mocks, maximum(D), binned_distances)
    if dt < d0
        return dt
    else
        mocks[random_point] = current
        return d0
    end
end

plot(progression, c=:black, lab="", dpi=400, lw=2)

plot(h_intra ./ sum(h_intra), m=:circle, c=:lightgrey, lab="Empirical distribution")
scatter!(h_faux ./ sum(h_faux), lab="Simulated data")