using Revise
using SimpleSDMLayers
using GBIF
using Plots
using Distances
using LinearAlgebra
using StatsBase
using KernelDensity

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
occurrences!(observations)

# Clip raster
elevation = clip(elevation, observations)

# Coordinates
coordinates = [(observations[i].longitude, observations[i].latitude) for i in 1:length(observations)]
_R = 6371.0
_dist = Distances.Haversine(_R)
intrasp_distances = pairwise(_dist, coordinates)
max_distance = maximum(vec(intrasp_distances))
h_intra = fit(Histogram, vec(intrasp_distances) ./ max_distance, 0.0:0.05:1.0).weights

# Cells within d
function randompoint(ref, d, R)
    λ, φ = deg2rad.(ref)
    δ = d / R

    α = deg2rad(rand() * 360.0)
    φ2 = asin(sin(φ) * cos(δ) + cos(φ) * sin(δ) * cos(α))
    λ2 = λ + atan(sin(α) * sin(δ) * cos(φ), cos(δ) - sin(φ) * sin(φ2))
    return rad2deg.((λ2, φ2))
end

# Initial point
invalid = true
while invalid
    # Get a random distance
    rdist = rand(intrasp_distances)
    rcent = rand(fauxpoints)
    new_faux = randompoint(rcent, rdist, _R)
    invalid = isnothing(elevation[new_faux...])
end
fauxpoints = [new_faux]

while length(fauxpoints) < length(coordinates)
    invalid = true
    while invalid
        # Get a random distance
        rdist = rand(intrasp_distances)
        rcent = rand(fauxpoints)
        new_faux = randompoint(rcent, rdist, _R)
        invalid = isnothing(elevation[new_faux...])
    end

    push!(fauxpoints, new_faux)
    fd = pairwise(_dist, fauxpoints)

    if maximum(vec(fd)) > max_distance
        @info "Point rejected (distance too high)"
        pop!(fauxpoints)
    end
end

progression = Float64[]
counter = 0
d0 = Inf

for i in 1:500_000
    rdist = rand(intrasp_distances)
    ridx = rand(1:length(fauxpoints))
    invalid = true
    while invalid
        # Get a random distance
        rdist = rand(intrasp_distances)
        rcent = fauxpoints[ridx]
        new_faux = randompoint(rcent, rdist, _R)
        invalid = isnothing(elevation[new_faux...])
    end
    oldpoint = fauxpoints[ridx]
    fauxpoints[ridx] = new_faux
    fd = pairwise(_dist, fauxpoints)
    h_faux = fit(Histogram, vec(fd) ./ max_distance, 0.0:0.05:1.0).weights
    dt = kl_divergence(h_faux ./ sum(h_faux), h_intra ./ sum(h_intra))
    if dt < d0
        d0 = dt
        @info "Score: $(d0)"
        counter = 0
    else
        fauxpoints[ridx] = oldpoint
        counter = counter + 1
    end
    push!(progression, d0)
    if counter == 100
        @info "No improvement for the last 100 steps, returning"
        break
    end
end

plot(progression, c=:black, lab="", dpi=400, lw=2)

plot(h_intra ./ sum(h_intra), m=:circle, c=:lightgrey, lab="Empirical distribution")
scatter!(h_faux ./ sum(h_faux), lab="Simulated data")