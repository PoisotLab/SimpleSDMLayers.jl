using Revise
using SimpleSDMLayers
using GBIF
using Plots
using Distances
using LinearAlgebra
using StatsBase
using BenchmarkTools
using ProgressMeter

Df = Distances.Haversine(6371.0)

"""
Get the coordinates for a list of observations, filtering the ones that do not
correspond to valid layer positions
"""
function coordinates(observations, layer)
    xy = [(observations[i].longitude, observations[i].latitude) for i in 1:length(observations)]
    filter!(c -> !isnothing(layer[c...]), xy)
    return hcat(collect.(unique(xy))...)
end

"""
Bin a distance matrix, where m is the maximum distance allowed
"""
function bin_distances(D, m)
    w = fit(Histogram, vec(D) ./ m, LinRange(0.0, 1.0, 21)).weights
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
layer - for ref
xy - points
Dxy - distances
"""
function new_random_point(layer, xy, Dxy)
    invalid = true
    global point
    while invalid
        point = randompoint(xy[:, rand(1:size(xy, 2))], rand(Dxy))
        invalid = isnothing(layer[point...])
    end
    return point
end

"""
Generates the initial proposition for points
"""
function generate_initial_points(layer, xy, Dxy)
    all_points = copy(xy)
    all_points[:, 1] .= new_random_point(layer, xy, Dxy)
    for i in 2:size(xy, 2)
        global point
        invalid = true
        while invalid
            point = new_random_point(layer, all_points[:, 1:(i-1)], Dxy)
            all_points[:, i] .= point
            invalid = maximum(pairwise(Df, all_points[:, 1:i])) > maximum(Dxy)
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
function distribution_distance(x, y)
    m = max(maximum(x), maximum(y))
    p = bin_distances(x, m)
    q = bin_distances(y, m)
    return sqrt(js_divergence(p, q) / log(2))
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

function fauxcurrence(layer, xy::Vector{Tuple{Float64,Float64}}; stop_at=0.01, max_iter=10_000)
    Dx = distance_matrix(xy)
    mocks = _initial_proposition(layer, xy, Dx)
    Dy = distance_matrix(mocks)
    d0 = _points_distance(Dx, Dy)
    for i in 1:max_iter
        d0 = _improve_one_point!(mocks, layer, Dx, d0)
        @info "Current d₀:\t$(progression[i])"
        if progression[i] < stop_at
            @info "Breaking after $(i) iterations"
            break
        end
    end
    return mocks
end


#=
# Hawaii example
_bbox = (left=-160.0, right=-154.5, bottom=18.5, top=22.5)
layer = convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; _bbox..., resolution=0.5))
plot(layer, frame=:box, c=:bamako, dpi=400)
taxa = [
    GBIF.taxon("Himatione sanguinea"; strict=true),
    GBIF.taxon("Paroaria capitata"; strict=true),
    GBIF.taxon("Pluvialis fulva"; strict=true),
    GBIF.taxon("Pandion haliaetus"; strict=true)
]
=#

_bbox = (left=-80.0, right=-56.0, bottom=44.0, top=62.0)
layer = convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; _bbox..., resolution=0.5))
plot(layer, frame=:box, c=:bamako, dpi=400)
taxa = [
    GBIF.taxon("Vulpes vulpes"; strict=true),
    GBIF.taxon("Urocyon cinereoargenteus"; strict=true),
    GBIF.taxon("Vulpes lagopus"; strict=true)
]

observations = []
for t in taxa
    obs = occurrences(t,
        "hasCoordinate" => "true",
        "decimalLatitude" => (_bbox.bottom, _bbox.top),
        "decimalLongitude" => (_bbox.left, _bbox.right),
        "limit" => 300)
    push!(observations, obs)
end

xy = [coordinates(obs, layer) for obs in observations]
Dxy = [pairwise(Df, xy[i]) for i in 1:length(xy)]
fc = [generate_initial_points(layer, xy[i], Dxy[i]) for i in 1:length(xy)]
Dfc = [pairwise(Df, xy[i], fc[i]) for i in 1:length(fc)]
JS = [distribution_distance(Dxy[i], Dfc[i]) for i in 1:length(fc)]
Dij = [pairwise(Df, xy[i], xy[j]) for i in 1:(length(xy)-1) for j in (i+1):length(xy)]
Dkl = [pairwise(Df, fc[i], fc[j]) for i in 1:(length(xy)-1) for j in (i+1):length(xy)]
JSp = [distribution_distance(Dij[i], Dkl[i]) for i in 1:length(Dij)]
optimum = mean(vcat(JS, JSp))

progress = zeros(Float64, 10_000)
progress[1] = optimum

for i in 2:length(progress)
    # Get a random set of points to change
    set_to_change = rand(1:length(fc))

    # Get a random point to change in the layer
    point_to_change = rand(1:size(fc[set_to_change], 2))

    # Save the old point
    current_point = fc[set_to_change][:, point_to_change]

    # Generate a new proposition
    fc[set_to_change][:, point_to_change] .= new_random_point(layer, current_point, Dxy[set_to_change])

    # Get the pairwise distance matrices
    Dfc = [pairwise(Df, xy[i], fc[i]) for i in 1:length(fc)]
    JS = [distribution_distance(Dxy[i], Dfc[i]) for i in 1:length(fc)]

    # Same with the pairwise inter-specific distances
    Dij = [pairwise(Df, xy[i], xy[j]) for i in 1:(length(xy)-1) for j in (i+1):length(xy)]
    Dkl = [pairwise(Df, fc[i], fc[j]) for i in 1:(length(xy)-1) for j in (i+1):length(xy)]
    JSp = [distribution_distance(Dij[i], Dkl[i]) for i in 1:length(Dij)]
    d0 = mean(vcat(JS, JSp))

    if d0 < optimum
        optimum = d0
        @info optimum
    else
        fc[set_to_change][:, point_to_change] .= current_point
    end
    progress[i] = optimum
end

plot(progress)

Dx = distance_matrix(xy)
Dy = distance_matrix(mocks)
m = max(maximum(Dx), maximum(Dy))
plot(bin_distances(Dx, m), dpi=600, lw=1.0, c=:white, lab="Measured", lc=:black, m=:circle)
scatter!(bin_distances(Dy, m), c=:black, lab="Simulated", ms=3, m=:diamond)
xaxis!("Distance bin", 1:20)
yaxis!("Density", (0, 0.5))

p = [plot(layer, frame=:grid, c=:grey, cbar=false, legend=:bottomleft, size=(700, 300), dpi=600) for i in 1:length(xy), j in 1:2]
c = distinguishable_colors(length(xy) + 4)[(end-length(xy)+1):end]
for i in 1:length(xy)
    scatter!(p[i, 1], xy[i][1, :], xy[i][2, :], lab="", ms=2, c=c[i], msw=0.0)
    scatter!(p[i, 2], fc[i][1, :], fc[i][2, :], lab="", ms=2, c=c[i], msw=0.0, m=:diamond)
end
plot(p..., layout=(2, length(xy)))
savefig("demo.png")
