module SSLTestMVStats
    using SimpleSDMLayers
    using MultivariateStats
    using NeutralLandscapes
    using Test

    layers = [SimpleSDMPredictor(rand(MidpointDisplacement(0.9), (50,50))) for i in 1:30]


    @assert typeof(transform(fit(PCA, layers),layers)) <: Vector{T} where T<:SimpleSDMLayer
    @assert typeof(transform(fit(PPCA, layers),layers)) <: Vector{T} where T<:SimpleSDMLayer
    # @assert typeof(transform(fit(KernelPCA, layers),layers)) <: Vector{T} where T<:SimpleSDMLayer

   # @assert typeof(transform(fit(Whitening, layers),layers)) <: Vector{T} where T<:SimpleSDMLayer

end
