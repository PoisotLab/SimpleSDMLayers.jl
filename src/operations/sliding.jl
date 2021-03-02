function haversine(lon1, lat1, lon2, lat2; R=6371.0)
    φ1 = lat1 * π/180.0
    φ2 = lat2 * π/180.0
    Δφ = (lat2-lat1) * π/180.0
    Δλ = (lon2-lon1) * π/180.0
    a = sin(Δφ/2.0)^2.0 + cos(φ1)*cos(φ2) * sin(Δλ)^2.0
    c = 2.0 * atan(sqrt(a), sqrt(1.0-a))
    return R*c
end

haversine(p1, p2; R=6371.0) = haversine(p1..., p2...; R=R)

"""
    slidingwindow(L::LT, f::FT, d::IT) where {LT <: SimpleSDMLayer, FT <: Function, IT <: Number}

This function will replace the value at any cell by applying the function `f` to
the array of cells values that are within a distance `d` (in kilometers) from
the focal cell. This is, for example, useful to use an average to smooth out the
layers. The distance is estimated using the haversine distance, assuming that
the radius of the Earth is 6371.0 km. This means that the size of the window
will vary a little bit across latitudes, but this is far better than using a
number of cells, which would have dramatic consequences near the poles.

It *always* returns a `SimpleSDMResponse`, and the cells containing `nothing`
will also not contain a value in the output. This is *different* from the
behavior of `coarsen`, which tends to expand the area of the layer in which we
have data.

This function is currently relatively slow. Performance improvements will arrive
at some point.
"""
function slidingwindow(layer::LT, f::FT, d::IT; threaded::Bool=true) where {LT <: SimpleSDMLayer, FT <: Function, IT <: Number}
    # We infer the return type from a call to the function on the first three elements
    return_type = typeof(f(collect(layer)[1:min(3, length(layer))]))

    # New layer using typed similar
    N = similar(layer, return_type)

    # Store latitudes and longitudes
    _lat, _lon = latitudes(layer), longitudes(layer)

    # Vector of all positions with a value
    filled_positions = CartesianIndices(layer.grid)[findall(!isnothing, layer.grid)]

    # We then filter in the occupied positions
    if threaded
        for pos in filled_positions
            neighbors = filter(p -> haversine((_lon[Tuple(pos)[2]], _lat[Tuple(pos)[1]]), (_lon[Tuple(p)[2]], _lat[Tuple(p)[1]])) < d, filled_positions)
            N.grid[pos] = f(layer.grid[neighbors])
        end
    else
        Threads.@threads for pos in filled_positions
            neighbors = filter(p -> haversine((_lon[Tuple(pos)[2]], _lat[Tuple(pos)[1]]), (_lon[Tuple(p)[2]], _lat[Tuple(p)[1]])) < d, filled_positions)
            N.grid[pos] = f(layer.grid[neighbors])
        end
    end

    # And we return the object
    return N
end
