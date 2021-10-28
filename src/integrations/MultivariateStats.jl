# WARNING this file is only loaded if MultivariateStats.jl is also active
# This all happens thanks to the Requires.jl package

"""
	MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where T <: SimpleSDMLayer

	Overloads the `fit` function from `MultivariateStats.jl`.
"""
function MultivariateStats.fit(
    a::PT,
    layers::Vector{T},
    kwargs...,
) where {
    T<:SimpleSDMLayer,
    PT<:Union{Type{MultivariateStats.PCA},Type{MultivariateStats.PPCA}},
}
    input = _make_input(layers)
    pca = MultivariateStats.fit(a, input, kwargs...)
    return pca
end

"""
	MultivariateStats.transform(proj::P, layers::Vector{T}) where {P,T<:SimpleSDMLayer}

	Overload of the `transform` function from `MultivariateStats.jl`. 
"""
function MultivariateStats.transform(proj::P, layers::Vector{T}) where {P,T<:SimpleSDMLayer}
    input = _make_input(layers)
    numsites = length(layers)

    A = MultivariateStats.projection(proj)
    outdim = MultivariateStats.outdim(proj)
    transformedvals = zeros(outdim, numsites)

    for loc = 1:numsites
        transformedvals[:, loc] = A' * vec(input[:, loc])
    end

    newlayers = [zeros(Float32, size(layers[begin])) for l = 1:outdim]
    for loc = 1:numsites
        for layer = 1:outdim
            newlayers[layer][loc] = transformedvals[layer, loc]
        end
    end

    newlayers = SimpleSDMPredictor.(newlayers, boundingbox(layers[begin])...)
    return map(f -> rescale(f, (0, 1)), newlayers)
end


"""
	_make_input(layers)

	Creates a matrix for input to PCA. 
	The resulting matrix has columns for each shared location in 
	the vector of layers, where the value of each column is the 
"""
function _make_input(layers)
    common_keys = reduce(∩, keys.(layers))
    return hcat([layer[common_keys] for layer in layers]...)
end
