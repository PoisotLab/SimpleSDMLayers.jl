# Declares the relevant provided datasets
provides(::Type{CHELSA}, ::Type{BioClim}) = true

# Where to store them?
_rasterpath(::Type{CHELSA}) = "CHELSA"

# How to read these data?
_readfunction(::Type{CHELSA}, ::Type{<:LayerDataset}) = SimpleSDMLayers.geotiff

# Layer names - inherited from WorldClim
layernames(::Type{CHELSA}, ::Type{BioClim}) = layernames(WorldClim, BioClim)
