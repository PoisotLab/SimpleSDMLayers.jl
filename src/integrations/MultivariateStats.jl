# WARNING this file is only loaded if MultivariateStats.jl is also active
# This all happens thanks to the Requires.jl package

"""
	MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where T <: SimpleSDMLayer

	Overloads the `fit` function from `MultivariateStats.jl`.
"""
function MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where {T<:SimpleSDMLayer}
    _layers_are_compatible(layers) || return ArgumentError("layers are not compatible")
    input = _make_input(layers)
    proj = MultivariateStats.fit(a, input', kwargs...)
    return proj
end

"""
	MultivariateStats.transform(proj::P, layers::Vector{T}) where {P,T<:SimpleSDMLayer}

	Overload of the `transform` function from `MultivariateStats.jl`. Here `proj` is a 
	and output object from `MultivariateStats.fit` (see above). 
"""
function MultivariateStats.transform(proj::P, layers::Vector{T}) where {P,T<:SimpleSDMLayer}
    A = MultivariateStats.projection(proj)
    outdim = MultivariateStats.outdim(proj)    
    newlayers = [similar(layers[begin]) for i in 1:outdim]
    common_keys = reduce(∩, keys.(layers))

    for key in common_keys
        pcaproj = A' * vcat([layer[key] for layer in layers]...)
        for i in 1:outdim
            newlayers[i][key] = pcaproj[i]
        end
    end
    return map(f->rescale(f, (0,1)),newlayers)
end
