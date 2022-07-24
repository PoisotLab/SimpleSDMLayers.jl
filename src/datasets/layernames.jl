layernames(::Type{CHELSA}, ::Type{BioClim}) = layernames(WorldClim, BioClim)


function layernames(::Type{TerraClimate}, ::Type{PrimaryClimateVariable})
    return (
        "Maximum temperature",
        "Minimum temperature",
        "Vapor pressure",
        "Precipitation accumulation",
        "Downward surface shortwave radiation",
        "Wind speed"
    )
end

function layernames(::Type{TerraClimate}, ::Type{SecondaryClimateVariable})
    return (
        "Reference evapostranspiration",
        "Runoff",
        "Actual evapostranspiration",
        "Climate water deficit",
        "Soil moisture",
        "Snow water equivalent",
        "Palmer drought severity index",
        "Vapor pressure deficit"
    )
end
