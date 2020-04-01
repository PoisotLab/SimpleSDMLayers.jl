import Base: size
import Base: stride
import Base: eachindex
import Base: getindex
import Base: setindex!
import Base: similar
import Base: copy
import Base: eltype
import Base: convert

"""
    Base.convert(::Type{SimpleSDMResponse}, p::T) where {T <: SimpleSDMPredictor}

Returns a response with the same grid and bounding box as the predictor.
"""
function Base.convert(::Type{SimpleSDMResponse}, p::T) where {T <: SimpleSDMPredictor}
   return copy(SimpleSDMResponse(p.grid, p.left, p.right, p.bottom, p.top))
end

"""
    Base.convert(::Type{SimpleSDMPredictor}, p::T) where {T <: SimpleSDMResponse}

Returns a predictor with the same grid and bounding box as the response.
"""
function Base.convert(::Type{SimpleSDMPredictor}, p::T) where {T <: SimpleSDMResponse}
   return copy(SimpleSDMPredictor(p.grid, p.left, p.right, p.bottom, p.top))
end

"""
    Base.convert(::Type{Matrix}, p::T) where {T <: SimpleSDMLayer}

Returns the grid as an array.
"""
function Base.convert(::Type{Matrix}, p::T) where {T <: SimpleSDMLayer}
   return copy(p.grid)
end

"""
    Base.eltype(p::T) where {T <: SimpleSDMLayer}

Returns the type of the values stored in the grid.
"""
function Base.eltype(p::T) where {T <: SimpleSDMLayer}
   return eltype(p.grid)
end

"""
    Base.size(p::T) where {T <: SimpleSDMLayer}

Returns the size of the grid.
"""
function Base.size(p::T) where {T <: SimpleSDMLayer}
   return size(p.grid)
end

"""
    Base.size(p::T, i...) where {T <: SimpleSDMLayer}

Returns the size of the grid alongside a dimension.
"""

function Base.size(p::T, i...) where {T <: SimpleSDMLayer}
   return size(p.grid, i...)
end

"""
    Base.stride(p::T; dims::Union{Nothing,Integer}=nothing) where {T <: SimpleSDMLayer}

Returns the stride, *i.e.* the length, of cell dimensions, possibly alongside a
side of the grid.
"""
function Base.stride(p::T; dims::Union{Nothing,Integer}=nothing) where {T <: SimpleSDMLayer}
   lon_stride = (p.right-p.left)/size(p, 2)/2.0
   lat_stride = (p.top-p.bottom)/size(p, 1)/2.0
   dims == nothing && return (lon_stride, lat_stride)
   dims == 1 && return lon_stride
   dims == 2 && return lat_stride
end

"""
    Base.eachindex(p::T) where {T <: SimpleSDMLayer}

Returns the index of the grid.
"""
function Base.eachindex(p::T) where {T <: SimpleSDMLayer}
   return eachindex(p.grid)
end

"""
Extracts a  value from a layer by its grid position.
"""
function Base.getindex(p::T, i::Int64) where {T <: SimpleSDMLayer}
   return p.grid[i]
end

"""
    Base.getindex(p::T, i::R, j::R) where {T <: SimpleSDMLayer, R <: UnitRange}

Extracts a series of positions in a layer, and returns a layer corresponding to
the result. This is essentially a way to rapidly crop a layer to a given subset
of its extent. The `i` and `j` arguments are `UnitRange`s (of `Integer`).

The layer returned by this function will have the same type as the layer passed
as its argument, but this can be changed using `convert`.
"""
function Base.getindex(p::T, i::R, j::R) where {T <: SimpleSDMLayer, R <: UnitRange}
   return T(
            p.grid[i,j],
            minimum(longitudes(p)[j])-stride(p)[1],
            maximum(longitudes(p)[j])+stride(p)[1],
            minimum(latitudes(p)[i])-stride(p)[2],
            maximum(latitudes(p)[i])+stride(p)[2]
           )
end

"""
Given a layer and a latitude, returns NaN if the latitude is outside the
range, or the grid index containing this latitude if it is within range
"""
function match_latitude(p::T, l::K) where {T <: SimpleSDMLayer, K <: AbstractFloat}
   l > p.top && return NaN
   l < p.bottom && return NaN
   return findmin(abs.(l .- latitudes(p)))[2]
end

"""
Given a layer and a longitude, returns NaN if the longitude is outside the
range, or the grid index containing this longitude if it is within range
"""
function match_longitude(p::T, l::K) where {T <: SimpleSDMLayer, K <: AbstractFloat}
   l > p.right && return NaN
   l < p.left && return NaN
   return findmin(abs.(l .- longitudes(p)))[2]
end

"""
    Base.getindex(p::T, longitude::K, latitude::K) where {T <: SimpleSDMLayer, K <: AbstractFloat}

Extracts the value of a layer at a given latitude and longitude. If values
outside the range are requested, will return `NaN`.
"""
function Base.getindex(p::T, longitude::K, latitude::K) where {T <: SimpleSDMLayer, K <: AbstractFloat}
   i = match_longitude(p, longitude)
   j = match_latitude(p, latitude)
   isnan(i) && return NaN
   isnan(j) && return NaN
   return p.grid[j, i]
end

"""
    Base.getindex(p::T, longitudes::Tuple{K,K}, latitudes::Tuple{K,K}) where {T <: SimpleSDMLayer, K <: AbstractFloat}

Extracts a series of positions in a layer, and returns a layer corresponding to
the result. This is essentially a way to rapidly crop a layer to a given subset
of its extent. The `longitudes` and `latitudes` arguments are tuples of floating
point values, representing the bounding box of the layer to extract.

The layer returned by this function will have the same type as the layer passed
as its argument.
"""
function Base.getindex(p::T, longitudes::Tuple{K,K}, latitudes::Tuple{K,K}) where {T <: SimpleSDMLayer, K <: AbstractFloat}
   imin, imax = [match_longitude(p, l) for l in longitudes]
   jmin, jmax = [match_latitude(p, l) for l in latitudes]
   any(isnan.([imin, imax, jmin, jmax])) && throw(ArgumentError("Unable to extract, coordinates outside of range"))
   return p[jmin:jmax, imin:imax]
end


"""
    Base.getindex(p1::T1, p2::T2) where {T1 <: SimpleSDMLayer, T2 <: SimpleSDMLayer}

Extract a layer based on a second layer. Note that the two layers must be
*compatible*, which is to say they must have the same bounding box and grid
size.
"""
function Base.getindex(p1::T1, p2::T2) where {T1 <: SimpleSDMLayer, T2 <: SimpleSDMLayer}
   SimpleSDMLayers.are_compatible(l1, l2)
   return p1[(p2.left, p2.right), (p2.bottom, p2.top)]
end

"""
     Base.setindex!(p::SimpleSDMResponse{T}, v::T, i...) where {T}

Changes the value of a cell, or a range of cells, as indicated by their grid
positions.
"""
function Base.setindex!(p::SimpleSDMResponse{T}, v::T, i...) where {T}
   @assert typeof(v) <: eltype(p.grid)
   p.grid[i...] = v
end

"""
    Base.setindex!(p::T, v, lon::Float64, lat::Float64) where {T <: SimpleSDMResponse}

Changes the values of the cell including the point at the requested latitude and
longitude.
"""
function Base.setindex!(p::SimpleSDMResponse{T}, v::T, lon::Float64, lat::Float64) where {T}
   i = match_longitude(p, lon)
   j = match_latitude(p, lat)
   p[j,i] = v
end

"""
    Base.similar(l::T) where {T <: SimpleSDMLayer}

Returns a `SimpleSDMResponse` of the same dimensions as the original layer, with
`NaN` in the same positions. The rest of the values are replaced by the output
of `zero(eltype(p.grid))`, which implies that there must be a way to get a zero
for the type. If not, the same result can always be achieved through the use of
`copy`, manual update, and `convert`.
"""
function Base.similar(l::T) where {T <: SimpleSDMLayer}
   emptygrid = similar(l.grid)
   for i in eachindex(emptygrid)
      emptygrid[i] = isnan(l.grid[i]) ? NaN : zero(eltype(l.grid))
   end
   return SimpleSDMResponse(emptygrid, l.left, l.right, l.bottom, l.top)
end

"""
    Base.copy(l::T) where {T <: SimpleSDMLayer}

Returns a new copy of the layer, which has the same type.
"""
function Base.copy(l::T) where {T <: SimpleSDMLayer}
   copygrid = copy(l.grid)
   return T(copygrid, copy(l.left), copy(l.right), copy(l.bottom), copy(l.top))
end
