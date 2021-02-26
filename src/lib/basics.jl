"""
    latitudes(layer::T) where {T <: SimpleSDMLayer}

Returns an iterator with the latitudes of the SDM layer passed as its argument.
This returns the latitude at the center of each cell in the grid.
"""
function latitudes(layer::T) where {T <: SimpleSDMLayer}
    grid_size = stride(layer, 2)
    return range(layer.bottom+grid_size; stop=layer.top-grid_size, length=size(layer, 1))
end

"""
    longitudes(layer::T) where {T <: SimpleSDMLayer}

Returns an iterator with the longitudes of the SDM layer passed as its argument.
This returns the longitudes at the center of each cell in the grid.
"""
function longitudes(layer::T) where {T <: SimpleSDMLayer}
    grid_size = stride(layer, 1)
    return range(layer.left+grid_size; stop=layer.right-grid_size, length=size(layer, 2))
end

"""
    _layers_are_compatible(l1::X, l2::Y) where {X <: SimpleSDMLayer, Y <: SimpleSDMLayer}

    Internal function to verify if layers are compatible, i.e. have the same
    size and bounding coordinates.

"""
function _layers_are_compatible(l1::X, l2::Y) where {X <: SimpleSDMLayer, Y <: SimpleSDMLayer}
    size(l1) == size(l2) || throw(ArgumentError("The layers have different sizes"))
    l1.top == l2.top || throw(ArgumentError("The layers have different top coordinates"))
    l1.left == l2.left || throw(ArgumentError("The layers have different left coordinates"))
    l1.bottom == l2.bottom || throw(ArgumentError("The layers have different bottom coordinates"))
    l1.right == l2.right || throw(ArgumentError("The layers have different right coordinates"))
end

function _layers_are_compatible(layers::Array{T}) where {T <: SimpleSDMLayer}
    all(layer -> _layers_are_compatible(layer, layers[1]), layers)
end
