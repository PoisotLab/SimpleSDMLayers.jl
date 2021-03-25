function _get_raster(::Type{CHELSA}, ::Type{BioClim}, layer::Integer)
    1 ≤ layer ≤ 19 || throw(ArgumentError("The layer must be between 1 and 19"))

    path = joinpath(SimpleSDMLayers._layers_assets_path, _rasterpath(CHELSA), _rasterpath(BioClim))
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

function _get_raster(::Type{CHELSA}, ::Type{BioClim}, mod::CMIP5, fut::RepresentativeConcentrationPathway, layer::Integer, year="2041-2060")
    1 ≤ layer ≤ 19 || throw(ArgumentError("The layer must be between 1 and 19"))

    path = joinpath(SimpleSDMLayers._layers_assets_path, _rasterpath(CHELSA), _rasterpath(BioClim), _rasterpath(mod), _rasterpath(fut), year)
    isdir(path) || mkpath(path)

    root = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/cmip5/$(year)/bio/"
    filename = "CHELSA_bio_mon_$(_rasterpath(mod))_$(_rasterpath(fut))_r1i1p1_g025.nc_$(layer)_$(year)_V1.2.tif"
    
    filepath = joinpath(path, filename)
    if !(isfile(filepath))
        Downloads.download(root * filename, filepath)
    end

    return filepath
end