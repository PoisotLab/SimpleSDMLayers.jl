using Revise
using SimpleSDMLayers
using GBIF
using Plots
using Distances
using LinearAlgebra
using StatsBase
using BenchmarkTools

# One layer
_bbox = (left=117.0, right=126.2, bottom=-6.1, top=3.3)
elevation = convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; _bbox..., resolution=0.5))
plot(elevation, frame=:box, c=:bamako, dpi=400)

# Occurrences
observations = occurrences(
    GBIF.taxon("Draco beccarii"; strict=true),
    "hasCoordinate" => "true",
    "decimalLatitude" => (_bbox.bottom, _bbox.top),
    "decimalLongitude" => (_bbox.left, _bbox.right),
    "limit" => 300,
)
occurrences!(observations)

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
    return fit(Histogram, vec(D) ./ m, LinRange(0.0, 1.0, 10)).weights .+ 1
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
            proposition = randompoint(rand(xy), rand(D))
            invalid = isnothing(layer[proposition...])
        end
        all_points[i] = proposition
    end
    return all_points
end

function _points_distance(mocks, m, b)
    mD = distance_matrix(mocks)
    mb = bin_distances(mD, m)
    return kl_divergence(mb, b)
end

function _improve_one_point!(mocks, layer, D, b, d0)
    random_point = rand(eachindex(mocks))
    current = mocks[random_point]
    global proposition
    invalid = true
    counter = 0
    while invalid
        counter += 1
        random_distance = rand(D)
        proposition = randompoint(mocks[random_point], random_distance)
        mocks[random_point] = proposition
        invalid = isnothing(layer[proposition...])
        if counter >= 10
            mocks[random_point] = current
            return d0
        end
    end
    dt = _points_distance(mocks, maximum(D), b)
    if dt < d0
        return dt
    else
        mocks[random_point] = current
        return d0
    end
end

layer = copy(elevation)
xy = coordinates(observations, layer)
D = distance_matrix(xy)
b = bin_distances(D, maximum(D))
m = maximum(D)
mocks = _initial_proposition(layer, xy, D)
d0 = _points_distance(mocks, m, b)

progression = zeros(Float64, 10_000)
progression[1] = d0
for i in 2:length(progression)
    progression[i] = _improve_one_point!(mocks, layer, D, b, progression[i-1])
    @info "Time $(i)\t$(progression[i])"
end

plot(progression, c=:black, lab="", dpi=400, lw=2)

plot(b, dpi=400, lw=0.0, fill=(0, :grey, 0.2), lab="Measured")
plot!(bin_distances(distance_matrix(mocks), m), c=:black, lab="Simulated")

plot(layer, frame=:box, dpi=600, c=:grey, cbar=false)
scatter!(xy, c=:black, lab="Occurrences")
scatter!(mocks, c=:white, lab="Fauxcurrences", alpha=0.7, m=:square, ms=2)
savefig("demo.png")
