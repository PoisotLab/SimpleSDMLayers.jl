# Declares the relevant provided datasets
provides(::Type{EarthEnv}, ::Type{LandCover}) = true
provides(::Type{EarthEnv}, ::Type{Topography}) = true
provides(::Type{EarthEnv}, ::Type{HabitatHeterogeneity}) = true

# Where to store them?
_rasterpath(::Type{EarthEnv}) = "EarthEnv"
_rasterpath(::Type{EarthEnv}, ::Type{LandCover}) = joinpath(_rasterpath(WorldClim), "LandCover")
_rasterpath(::Type{EarthEnv}, ::Type{HabitatHeterogeneity}) = joinpath(_rasterpath(WorldClim), "HabitatHeterogeneity")
_rasterpath(::Type{EarthEnv}, ::Type{Topography}) = joinpath(_rasterpath(WorldClim), "Topography")

# Layer names
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