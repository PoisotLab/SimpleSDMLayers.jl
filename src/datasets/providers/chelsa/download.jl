function _get_raster(::Type{CHELSA}, ::Type{BioClim}, layer::Integer)
    1 ≤ layer ≤ 19 || throw(ArgumentError("The layer must be between 1 and 19"))

    path = joinpath(SimpleSDMLayers._layers_assets_path, _rasterpath(CHELSA), _rasterpath(BioClim))
    isdir(path) || mkpath(path)

    url_root = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/"
    filename = "CHELSA_bio$(layer)_1981-2010_V.2.1.tif"

    filepath = joinpath(path, filename)
    if !(isfile(filepath))
        Downloads.download(url_root * filename, filepath)
    end
    return filepath
end

function _get_raster(::Type{CHELSA}, ::Type{BioClim}, mod::CMIP6, fut::SharedSocioeconomicPathway, layer::Integer, year="2011-2040")
    @assert mod in [GFDLESM4, IPSLCM6ALR, MPIESM12HR, MRIESM20, UKESM10LL]
    1 ≤ layer ≤ 19 || throw(ArgumentError("The layer must be between 1 and 19"))
    year in ["2011-2040", "2041-2070", "2071-2100"] || throw(ArgumentError("The year must be 2011-2040, 2041-2070, or 2071-2100"))

    path = joinpath(SimpleSDMLayers._layers_assets_path, _rasterpath(CHELSA), _rasterpath(BioClim), _rasterpath(mod), _rasterpath(fut), year)
    isdir(path) || mkpath(path)

    root = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/$(year)/$(_rasterpath(mod))/$(_rasterpath(fut))/bio/"
    filename = "CHELSA_bio$(layer)_$(year)_$(lowercase(_rasterpath(mod)))_$(_rasterpath(fut))_V.2.1.tif"
    
    filepath = joinpath(path, filename)
    if !(isfile(filepath))
        Downloads.download(root * filename, filepath)
    end

    return filepath
end

