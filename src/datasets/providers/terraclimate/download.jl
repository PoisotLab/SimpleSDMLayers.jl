function _get_raster(::Type{TerraClimate}, ::Type{PrimaryClimateVariable}, layer::Integer; year=1958)
    _layer_max = length(layernames(TerraClimate, PrimaryClimateVariable))
    1 ≤ layer ≤ _layer_max || throw(ArgumentError("The layer must be between 1 and $(_layer_max)"))
    1958 ≤ year ≤ 2021 || throw(ArgumentError("The year must be between 1958 and 2021"))

    path = joinpath(SimpleSDMLayers._layers_assets_path, _rasterpath(TerraClimate), _rasterpath(PrimaryClimateVariable))
    isdir(path) || mkpath(path)

    # NOTE these are ordered as in layernames and this is very fragile
    layercodes = [
        "tmax", "tmin", "vp", "ppt", "srad", "ws"
    ]


    filename = "$(layercodes[layer])_$(year).nc"
    url_root = "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_"

    filepath = joinpath(path, filename)
    if !(isfile(filepath))
        Downloads.download(url_root * filename, filepath)
    end
    return filepath
end