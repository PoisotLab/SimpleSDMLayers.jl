layernames(::Type{CHELSA}, ::Type{BioClim}) = layernames(WorldClim, BioClim)

function layernames(::Type{EarthEnv}, ::Type{HabitatHeterogeneity})
    return (
        "Coefficient of variation",
        "Evenness",
        "Range",
        "Shannon",
        "Simpson",
        "Standard deviation",
        "Contrast",
        "Correlation",
        "Dissimilarity",
        "Entropy",
        "Homogeneity",
        "Maximum",
        "Uniformity",
        "Variance",
    )
end

function layernames(::Type{EarthEnv}, ::Type{LandCover})
    return (
        "Evergreen/Deciduous Needleleaf Trees",
        "Evergreen Broadleaf Trees",
        "Deciduous Broadleaf Trees",
        "Mixed/Other Trees",
        "Shrubs",
        "Herbaceous Vegetation",
        "Cultivated and Managed Vegetation",
        "Regularly Flooded Vegetation",
        "Urban/Built-up",
        "Snow/Ice",
        "Barren",
        "Open Water",
    )
end

function layernames(::Type{EarthEnv}, ::Type{Topography})
    return (
        "Elevation",
        "Slope",
        "Aspect Cosine",
        "Aspect Sine",
        "Aspect Eastness",
        "Aspect Northness",
        "Roughness",
        "Topographic Position Index",
        "Terrain Ruggedness Index",
        "Vector Ruggedness Measure",
        "∂(E-W slope)",
        "∂²(E-W slope)",
        "∂(N-S slope)",
        "∂²(N-S slope)",
        "Profile curvature",
        "Tangential curvature"
    )
end

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
