# WARNING this file is only loaded if MultivariateStats.jl is also active
# This all happens thanks to the Requires.jl package

"""
	MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where T <: SimpleSDMLayer

	Overloads the `fit` function from `MultivariateStats.jl`.
"""
function MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where {T<:SimpleSDMLayer}
    _layers_are_compatible(layers) || return ArgumentError("layers are not compatible")
    common_keys = reduce(∩, keys.(layers))
    input = hcat([vcat([layer[key] for layer in layers]...) for key in common_keys]...)
    proj = MultivariateStats.fit(a, input, kwargs...)
    return proj
end

"""
    transform(proj, layers::Vector{V}, 
    kwargs...) where {PT<:Union{MultivariateStats.PCA, MultivariateStats.PPCA},V<:SimpleSDMLayer} 

	Overload of the `transform` function from `MultivariateStats.jl`. Here `proj` is a 
	and output object from `MultivariateStats.fit` (see above). 
"""
function transform(
    proj,
    layers::AbstractVecOrMat{U},
    kwargs...,
) where {T<:Union{MultivariateStats.Whitening,MultivariateStats.PCA},U<:SimpleSDMLayer}
    outdim = MultivariateStats.outdim(proj)
    newlayers = [similar(layers[begin]) for i = 1:outdim]
    common_keys = reduce(∩, keys.(layers))

    input = hcat([vcat([layer[key] for layer in layers]...) for key in common_keys]...)

    for (ct, key) in enumerate(common_keys)
        pcaproj = MultivariateStats.transform(proj, input[:, ct])
        for i = 1:outdim
            newlayers[i][key] = pcaproj[i]
        end
    end
    return newlayers
end
