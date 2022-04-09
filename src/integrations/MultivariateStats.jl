# WARNING this file is only loaded if MultivariateStats.jl is also active
# This all happens thanks to the Requires.jl package

_allowed_transforms = (MultivariateStats.PCA, MultivariateStats.PPCA, MultivariateStats.KernelPCA, MultivariateStats.Whitening)
AllowedMultivariateTransforms = Union{_allowed_transforms...}

"""
	MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where T <: SimpleSDMLayer

Overloads the `fit` function from `MultivariateStats.jl`.
"""
function MultivariateStats.fit(proj::Type{K}, layers::Vector{T}; kwargs...) where {T<:SimpleSDMLayer, K <: AllowedMultivariateTransforms}
    _layers_are_compatible(layers) || return ArgumentError("layers are not compatible")
    common_keys = reduce(∩, keys.(layers))
    input = hcat([vcat([layer[key] for layer in layers]...) for key in common_keys]...)
    proj = MultivariateStats.fit(proj, input; kwargs...)
    return proj
end

"""
    transform(proj, layers::Vector{V}, 
    kwargs...) where {PT<:Union{MultivariateStats.PCA, MultivariateStats.PPCA},V<:SimpleSDMLayer} 

Overload of the `transform` function from `MultivariateStats.jl`. Here `proj` is an output object from `MultivariateStats.fit` (see above). 
"""
function MultivariateStats.transform(
    proj, layers::AbstractVecOrMat{U}; kwargs...
) where {U<:SimpleSDMLayer}
    _layers_are_compatible(layers) || return ArgumentError("layers are not compatible")
    
    newlayers = [similar(first(layers)) for i in 1:MultivariateStats.outdim(proj)]
    common_keys = reduce(∩, keys.(layers))

    input = hcat([vcat([layer[key] for layer in layers]...) for key in common_keys]...)

    for (ct, key) in enumerate(common_keys)
        pcaproj = MultivariateStats.transform(proj, input[:, ct])
        for i in 1:length(newlayers)
            newlayers[i][key] = pcaproj[i]
        end
    end

    return newlayers
end
