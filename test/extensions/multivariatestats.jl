module SSLTestMVStats
using SimpleSDMLayers
using MultivariateStats
using NeutralLandscapes
using Test

TEST_DIMS = (150, 150)
TEST_AUTOCORRELATION = 0.9
TEST_NUM_LAYERS = 10

layers = [
    SimpleSDMResponse(rand(MidpointDisplacement(TEST_AUTOCORRELATION), TEST_DIMS...)) for
    i in 1:TEST_NUM_LAYERS
]

@test typeof(transform(fit(PCA, layers), layers)) <: Vector{T} where {T<:SimpleSDMLayer}
#@test typeof(transform(fit(PPCA, layers),layers)) <: Vector{T} where T<:SimpleSDMLayer

#@assert typeof(transform(fit(KernelPCA, layers),layers)) <: Vector{T} where T<:SimpleSDMLayer
#@assert typeof(transform(fit(Whitening, layers),layers)) <: Vector{T} where T<:SimpleSDMLayer

end
