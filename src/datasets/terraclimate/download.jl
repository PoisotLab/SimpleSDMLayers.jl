function _get_raster(::Type{TerraClimate}, ::Type{PrimaryClimateVariable}, layer::Integer)
    1 ≤ layer ≤ length(layernames(TerraClimate, PrimaryClimateVariable)) || throw(ArgumentError("The layer must be between 1 and 19"))

    path = joinpath(SimpleSDMLayers._layers_assets_path, _rasterpath(TerraClimate), _rasterpath(PrimaryClimateVariable))
    isdir(path) || mkpath(path)
    
    layer = lpad(layer, 2, "0")
    filename = "CHELSA_bio10_$(layer).tif"
    url_root = "https://envicloud.os.zhdk.cloud.switch.ch/chelsa/chelsa_V1/climatologies/bio/"

    filepath = joinpath(path, filename)
    if !(isfile(filepath))
        Downloads.download(url_root * filename, filepath)
    end
    return filepath
end