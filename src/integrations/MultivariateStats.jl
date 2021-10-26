# WARNING this file is only loaded if MultivariateStats.jl is also active
# This all happens thanks to the Requires.jl package

"""
	MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where T <: SimpleSDMLayer
"""
function MultivariateStats.fit(a, layers::Vector{T}, kwargs...) where T <: SimpleSDMLayer
	input = _make_input(layers)
	pca = MultivariateStats.fit(a, input, kwargs...)
	return pca
end

"""
	#
"""
function MultivariateStats.transform(proj::P, layers::Vector{T}) where {P,T<:SimpleSDMLayer}
	input = _make_input(layers)
	numsites = length(layers)

	A = MultivariateStats.projection(proj)
	outdim = MultivariateStats.outdim(proj)
	transformedvals = zeros(outdim, numsites)
	
	for loc in 1:numsites
		transformedvals[:,loc] =  A' * vec(input[:,loc])
	end

	newlayers = [zeros(Float32, size(layers[begin])) for l in 1:outdim]
	for loc in 1:numsites
		for layer in 1:outdim
			newlayers[layer][loc] =  transformedvals[layer,loc]
		end
	end

	newlayers = SimpleSDMPredictor.(newlayers, boundingbox(layers[begin])...)
	newlayers = map(f -> rescale(f, (0,1)), newlayers)

	newlayers
end

function _make_input(layers)
	common_keys = reduce(∩, keys.(layers))
	x = hcat([layer[common_keys] for layer in layers]...);
	x
end

