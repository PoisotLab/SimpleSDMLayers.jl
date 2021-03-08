"""
    mosaic(f::TF, layers::Vector{T}, ::RT; sanitize::Bool=true, rtype::Type=Float64) where {TF <: Function, T <: SimpleSDMLayer}

Joins a series of _possibly_ overlapping `layers` by applying the function `f`
to the values that occupy the same cells. Note that the function `f`  should
return a single value and accept an vector as input. Functions like
`Statistics.mean`, etc, work well. `RT` is the return type of the function `f`.

The `sanitize` keyword is used to determine whether values that are `NaN` or
`nothing` should be removed before passing the vector of values to `f`. Using
`sanitize=false` can be useful if you want to do things like finding which layer
has the maximal value, in a context where some layers may have `nothing`.
"""
function mosaic(f::TF, layers::Vector{T}, ::Type{RT}; sanitize::Bool=true) where {TF <: Function, T <: SimpleSDMLayer, RT}

    # Check the dimensions
    for i in 1:(length(layers)-1)
        for j in 1:length(layers)
            if !(all(stride(layers[i]) .≈ stride(layers[j])))
                throw(DimensionMismatch("Layers $i and $j have different strides"))
            end
        end
    end

    # Get the new bounding boxes
    n_left = minimum([layer.left for layer in layers])
    n_right = maximum([layer.right for layer in layers])
    n_bottom = minimum([layer.bottom for layer in layers])
    n_top = maximum([layer.top for layer in layers])
    
    # Get the gridsize
    nr = round(Int64, (n_top - n_bottom)/2stride(layers[1],1))
    nc = round(Int64, (n_right - n_left)/2stride(layers[1],2))
    
    # Prepare the grid
    grid = fill(nothing, nc, nr)
    grid = convert(Matrix{Union{Nothing,RT}}, grid)
    L = SimpleSDMResponse(grid, n_left, n_right, n_bottom, n_top)

    # Fill in the information
    for lat in latitudes(L)
        for lon in longitudes(L)
            V = Any[layer[lon, lat] for layer in layers]
            if sanitize
                filter!(!isnothing, V)
                filter!(!isnan, V)
            end
            length(V) == 0 && continue
            r = f(V)
            L[lon, lat] = isnothing(r) ? nothing : convert(RT, r)
        end
    end

    # Return
    return T <: SimpleSDMResponse ? L : convert(SimpleSDMPredictor, L)
end