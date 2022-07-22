# Declares the relevant provided datasets
provides(::Type{WorldClim}, ::Type{BioClim}) = true
provides(::Type{WorldClim}, ::Type{Elevation}) = true

# Names for the layers in each dataset
function layernames(::Type{WorldClim}, ::Type{BioClim})
    return (
        "Annual Mean Temperature",
        "Mean Diurnal Range",
        "Isothermality",
        "Temperature Seasonality",
        "Max Temperature of Warmest Month",
        "Min Temperature of Coldest Month",
        "Temperature Annual Range",
        "Mean Temperature of Wettest Quarter",
        "Mean Temperature of Driest Quarter",
        "Mean Temperature of Warmest Quarter",
        "Mean Temperature of Coldest Quarter",
        "Annual Precipitation",
        "Precipitation of Wettest Month",
        "Precipitation of Driest Month",
        "Precipitation Seasonality",
        "Precipitation of Wettest Quarter",
        "Precipitation of Driest Quarter",
        "Precipitation of Warmest Quarter",
        "Precipitation of Coldest Quarter",
    )
end

layernames(::Type{WorldClim}, ::Type{Elevation}) = ("Elevation",)