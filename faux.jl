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
function bin_distances(D::Matrix{Float64}, m::Float64)::Vector{Float64}
    bins = 20
    r = LinRange(0.0, m+eps(typeof(m)), bins+1)
    c = zeros(eltype(D), bins)
    for i in 1:bins
        c[i] = count(D .< r[i+1])
        if i > 1
            c[i] = c[i] - sum(c[1:(i-1)])
        end
    end
    return c./sum(c)
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

#=
_bbox = (left=-80.0, right=-65.0, bottom=44.0, top=50.0)
layer = convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; _bbox..., resolution=0.5))
plot(layer, frame=:box, c=:bamako, dpi=400)
taxa = [
    GBIF.taxon("Canis latrans"; strict=true),
    GBIF.taxon("Alces alces"; strict=true),
    GBIF.taxon("Didelphis virginiana"; strict=true)
]
=#

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

observations = []
for t in taxa
    obs = occurrences(t,
        "hasCoordinate" => "true",
        "decimalLatitude" => (_bbox.bottom, _bbox.top),
        "decimalLongitude" => (_bbox.left, _bbox.right),
        "limit" => 300)
    occurrences!(obs)
    occurrences!(obs)
    push!(observations, obs)
end

# Generate the observation distances
obs = [coordinates(obs, layer) for obs in observations]
obs_intra_matrices = [zeros(Float64, (size(obs[i], 2), size(obs[i], 2))) for i in 1:length(obs)]
obs_inter_matrices = [zeros(Float64, (size(obs[i], 2), size(obs[j], 2))) for i in 1:(length(obs)-1) for j in (i+1):length(obs)]
for i in 1:length(obs)
    pairwise!(obs_intra_matrices[i], Df, obs[i])
end
cursor = 1
for i in 1:(length(obs)-1)
    for j in (i+1):length(obs)
        pairwise!(obs_inter_matrices[cursor], Df, obs[i], obs[j])
        cursor += 1
    end
end

# Generate the simulation distances
sim = [generate_initial_points(layer, obs[i], obs_intra_matrices[i]) for i in 1:length(obs)]
sim_intra_matrices = [zeros(Float64, (size(obs[i], 2), size(obs[i], 2))) for i in 1:length(obs)]
sim_inter_matrices = [zeros(Float64, (size(obs[i], 2), size(obs[j], 2))) for i in 1:(length(obs)-1) for j in (i+1):length(obs)]
for i in 1:length(sim)
    pairwise!(sim_intra_matrices[i], Df, sim[i])
end
cursor = 1
for i in 1:(length(sim)-1)
    for j in (i+1):length(sim)
        pairwise!(sim_inter_matrices[cursor], Df, sim[i], sim[j])
        cursor += 1
    end
end

JSintra = [distribution_distance(obs_intra_matrices[i], sim_intra_matrices[i]) for i in 1:length(sim_intra_matrices)]
JSinter = [distribution_distance(obs_inter_matrices[i], sim_inter_matrices[i]) for i in 1:length(sim_inter_matrices)]
JS = vcat(JSintra, JSinter)
optimum = mean(JS)

progress = zeros(Float64, 10_000)
scores = zeros(Float64, (length(JS), length(progress)))
scores[:, 1] .= JS
progress[1] = optimum

@profview for i in 2:length(progress)
    # Get a random set of points to change
    set_to_change = rand(1:length(sim))

    # Get a random point to change in the layer
    point_to_change = rand(1:size(sim[set_to_change], 2))

    # Save the old point
    current_point = sim[set_to_change][:, point_to_change]

    # Generate a new proposition
    sim[set_to_change][:, point_to_change] .= new_random_point(layer, current_point, obs_intra_matrices[set_to_change])

    # Update the distance matrices
    pairwise!(sim_intra_matrices[set_to_change], Df, sim[set_to_change])
    cursor = 1
    for i in 1:(length(sim)-1)
        for j in (i+1):length(sim)
            # Only update the matrices for which there was a change
            if (i == set_to_change)|(j == set_to_change)
                pairwise!(sim_inter_matrices[cursor], Df, sim[i], sim[j])
            end
            cursor += 1
        end
    end

    JSintra = [distribution_distance(obs_intra_matrices[i], sim_intra_matrices[i]) for i in 1:length(sim_intra_matrices)]
    JSinter = [distribution_distance(obs_inter_matrices[i], sim_inter_matrices[i]) for i in 1:length(sim_inter_matrices)]
    JS = vcat(JSintra, JSinter)

    d0 = mean(JS)
    scores[:, i] .= JS

    if d0 < optimum
        optimum = d0
        @info "t = $(lpad(i, 8)) λ = $(round(optimum; digits=4))"
    else
        sim[set_to_change][:, point_to_change] .= current_point
    end
    progress[i] = optimum
end

# Performance change plot
plot(scores', lab="", frame=:box)
plot!(progress, c=:black, lw=2, lab="Overall")
xaxis!("Iteration step")
yaxis!("Absolute performance")

# Map
p = [plot(layer, frame=:grid, c=:grey, cbar=false, legend=:bottomleft, size=(700, 300), dpi=600) for i in 1:length(xy), j in 1:2]
c = distinguishable_colors(length(xy) + 4)[(end-length(xy)+1):end]
for i in 1:length(xy)
    scatter!(p[i, 1], xy[i][1, :], xy[i][2, :], lab="", ms=2, c=c[i], msw=0.0)
    scatter!(p[i, 2], fc[i][1, :], fc[i][2, :], lab="", ms=2, c=c[i], msw=0.0, m=:diamond)
end
plot(p..., layout=(2, length(xy)))
savefig("demo.png")
