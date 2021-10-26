# WARNING this file is only loaded if MultivariateStats.jl is also active
# This all happens thanks to the Requires.jl package

function MultivariateStats.fit(PCA, layers::Vector{T}) where T <: SimpleSDMLayer
	_make_pca_layers(layers)
end

function transform(W, layers::Vector{T}) where T<: SimpleSDMLayer end
function transform!(W, layers::Vector{T}) where T<: SimpleSDMLayer end

function _make_pca_input(layers)
	numdims = length(layers)
	numsites = prod(size(layers[begin]))
	pcainput = zeros(numdims, numsites)

	for i in eachindex(layers[1])
		for layernum in 1:length(layers)
			if !isnothing(layers[layernum][i])
				pcainput[layernum,i] = layers[layernum][i]
			else
			end
		end
	end
	pcainput
end

function _make_pca_layers(layers)
	pcainput = _make_pca_input(layers)
	numdims, numsites = length(layers), prod(size(layers[begin]))

	pca = MultivariateStats.fit(MultivariateStats.PCA, pcainput, pratio=0.995)
	A = MultivariateStats.projection(pca)

	outdim = MultivariateStats.outdim(pca)

	transformedvals = zeros(outdim, numsites)
	
	for loc in 1:numsites
		transformedvals[:,loc] =  A' * vec(pcainput[:,loc])
	end


	newlayers = [zeros(Float32, size(layers[begin])) for l in 1:outdim]

	for loc in 1:numsites
		for layer in 1:outdim
			newlayers[layer][loc] =  transformedvals[layer,loc]
		end
	end


	newlayers = SimpleSDMPredictor.(newlayers, boundingbox(layers[begin])...)
	newlayers = map(f -> rescale(f, (0,1)), newlayers)
end
