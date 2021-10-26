module SSLTestMVStats
    using SimpleSDMLayers
    using MultivariateStats
    using NeutralLandscapes
    using Test

    layers = [SimpleSDMPredictor(rand(MidpointDisplacement(0.5), (50,50))) for i in 1:30]

    @assert typeof(fit(PCA, layers)) <: Vector{T} where T<:SimpleSDMLayer

end
