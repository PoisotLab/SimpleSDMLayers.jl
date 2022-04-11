using Revise
using SimpleSDMLayers
using GBIF
using Plots
using Distances
using LinearAlgebra
using StatsBase
using BenchmarkTools

# One layer
_bbox = (left=-160., right=-154.5, bottom=18.5, top=22.5)
layer = convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; _bbox..., resolution=0.5))
plot(layer, frame=:box, c=:bamako, dpi=400)

# Occurrences
observations = occurrences(
    GBIF.taxon("Himatione sanguinea"; strict=true),
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
    w = fit(Histogram, vec(D)./m, LinRange(0.0, 1.0, 20)).weights
    return w ./ sum(w)
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
            all_points[i] = proposition
            invalid = isnothing(layer[proposition...]) || (maximum(distance_matrix(all_points[1:i])) > maximum(D))
        end
    end
    return all_points
end

"""
Returns the Jensen-Shannon distance (i.e. the square root of the divergence) for
the two distance matrices. This version is prefered to the KL divergence in the
original implementation as it prevents the `Inf` values when p(x)=0 and q(x)>0.
The JS divergences is bounded between 0 and the natural log of 2, which gives an
absolute measure of fit allowing to compare the solutions. Note that the value
returned is *already* corrected, so it can be at most 1.0, and at best
(identical matrices) 0.
"""
function _points_distance(x, y)
    m = max(maximum(x), maximum(y))
    p = bin_distances(x, m)
    q = bin_distances(y, m)
    return sqrt(js_divergence(p, q)/log(2))
end

function _improve_one_point!(mocks, layer, D, d0)
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
    dt = _points_distance(D, distance_matrix(mocks))
    if dt < d0
        return dt
    else
        mocks[random_point] = current
        return d0
    end
end

xy = coordinates(observations, layer)
Dx = distance_matrix(xy)
mocks = _initial_proposition(layer, xy, Dx)
Dy = distance_matrix(mocks)
d0 = _points_distance(Dx, Dy)

progression = zeros(Float64, 10_000)
progression[1] = d0
for i in 2:length(progression)
    progression[i] = _improve_one_point!(mocks, layer, Dx, progression[i-1])
    @info "Time $(i)\t$(progression[i])"
end

plot(progression, c=:black, lab="", dpi=400, lw=2, ylab="Absolute fit", xlab="Epoch")

Dx = distance_matrix(xy)
Dy = distance_matrix(mocks)
m = max(maximum(Dx), maximum(Dy))
plot(bin_distances(Dx, m), dpi=400, lw=0.0, fill=(0, :grey, 0.2), lab="Measured")
plot!(bin_distances(Dy, m), c=:black, lab="Simulated")

plot(layer, frame=:box, dpi=600, c=:grey, cbar=false, legend=:bottomleft)
scatter!(xy, c=:black, lab="Occurrences")
scatter!(mocks, c=:white, lab="Fauxcurrences", alpha=0.7, m=:square, ms=2)
savefig("demo.png")
