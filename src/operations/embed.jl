"""
    This file contains interfaces to MVStats for 
    taking many SDMLayers and embedding their values
    in a lower dimensional space to create a smaller 
    number of SDMLayers

"""

function embed(layers::Vector{T}, PCA) where T <: SimpleSDMLayer

end

function make_pca_input(layers) 
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

function make_pca_layers(layers)
	pcainput = make_pca_input(layers)
	numdims, numsites = length(layers), prod(size(layers[begin]))

	pca = fit(PCA, pcainput, pratio=0.995)
	A = projection(pca)
	
	transformedvals = zeros(outdim(pca), numsites)
	
	for loc in 1:numsites
		transformedvals[:,loc] =  A' * vec(pcainput[:,loc])
	end
	
	
	### TODO this has to be  a similar SDMlayer to retain coordinates
	newlayers = [zeros(Float32, size(layers[begin])) for l in 1:outdim(pca)]	
	
	for loc in 1:numsites
		for layer in 1:outdim(pca)
			newlayers[layer][loc] =  transformedvals[layer,loc]
		end
	end
	
	
	newlayers = SimpleSDMPredictor.(newlayers, boundingbox(layers[begin])...)
	newlayers = map(f -> rescale(f, (0,1)), newlayers)
end