
function netcdf(::Type{LT}, filename, variablename, layerid::Integer=1; left=nothing, right=nothing, bottom=nothing, top=nothing) where {LT<:SimpleSDMLayer}
    # The file must exist
    @assert isfile(filename)

    # The layer for time must be valid
    @assert 1 <= layerid <= length(NetCDF.ncread(filename, "time"))

    # Coordinates for the part we want to extract
    lon_stride = NetCDF.ncgetatt(filename, "global", "geospatial_lon_resolution")
    lat_stride = NetCDF.ncgetatt(filename, "global", "geospatial_lat_resolution")
    lon_min = NetCDF.ncgetatt(filename, "global", "geospatial_lon_min")
    lon_max = NetCDF.ncgetatt(filename, "global", "geospatial_lon_max")
    lat_min = NetCDF.ncgetatt(filename, "global", "geospatial_lat_min")
    lat_max = NetCDF.ncgetatt(filename, "global", "geospatial_lat_max")

    # Correct the left/right/bottom/top values
    left = isnothing(left) ? lon_min - 0.5lon_stride : max(left, lon_min - 0.5lon_stride)
    right = isnothing(right) ? lon_max + 0.5lon_stride : min(right, lon_max + 0.5lon_stride)
    bottom = isnothing(bottom) ? lat_min - 0.5lat_stride : max(bottom, lat_min - 0.5lat_stride)
    top = isnothing(top) ? lat_max + 0.5lat_stride : min(top, lat_max + 0.5lat_stride)

    # Get the array positions to read
    left_idx = findfirst(l -> left <= l + 0.5lon_stride, NetCDF.ncread(filename, "lon"))
    right_idx = findlast(l -> right >= l - 0.5lon_stride, NetCDF.ncread(filename, "lon"))
    # BEWARE the Northernmost latitudes are at the top!
    bottom_idx = findlast(l -> bottom <= l - 0.5lat_stride, NetCDF.ncread(filename, "lat"))
    top_idx = findfirst(l -> top >= l + 0.5lat_stride, NetCDF.ncread(filename, "lat"))

    # Extract and rescale the values
    content = NetCDF.ncread(filename, variablename, start=[left_idx, top_idx, layerid], count=[right_idx - left_idx + 1, bottom_idx - top_idx + 1, 1])[:, :, 1]
    scale_factor = NetCDF.ncgetatt(filename, variablename, "scale_factor")
    content = content .* scale_factor

    # Get the missing value
    missing_value = NetCDF.ncgetatt(filename, variablename, "missing_value")

    # Figure out the correct bounding box
    _bbox = (
        left=NetCDF.ncread(filename, "lon")[left_idx] - 0.5lon_stride,
        right=NetCDF.ncread(filename, "lon")[right_idx] + 0.5lon_stride,
        bottom=NetCDF.ncread(filename, "lat")[bottom_idx] - 0.5lat_stride,
        top=NetCDF.ncread(filename, "lat")[top_idx] + 0.5lat_stride
    )

    # Prepare the layer
    content = convert(Matrix{Union{Nothing,eltype(content)}}, content)
    replace!(content, missing_value * scale_factor => nothing)

    return LT(rotl90(content), _bbox...)

end