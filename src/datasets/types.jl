"""
    WorldClim

[WorldClim](https://worldclim.org/) offers bioclimatic data both historical and
future, under `CMIP6` scenarios.

This provider currently offers `BioClim` data, both historical and future under
`CMIP6`, and the `Elevation` raster.
"""
struct WorldClim <: LayerProvider end

"""
    CHELSA

[CHELSA](https://chelsa-climate.org/) offers high resolution climatologies. This
provider currently offers `BioClim` data, both historical and future under
`CMIP5` *and* `CMIP6`.

Note that CHELSA offers a subset of all possible CMIP6 combinations, which is
supposed to be the most informative.
"""
struct CHELSA <: LayerProvider end

"""
    EarthEnv

Data from the [EarthEnv](https://www.earthenv.org//) project, all released under
a CC-BY-NC license to Tuanmu & Jetz. This provider currently offers `LandCover`
and `HabitatHeterogeneity` rasters.
"""
struct EarthEnv <: LayerProvider end

"""
    TerraClimate

Data from `https://www.climatologylab.org/terraclimate.html`, offering the
`PrimaryClimate` and `DerivedClimate` layers. The `TerraClimate` dataset offers
monthly variables.
"""
struct TerraClimate <: LayerProvider end

struct PrimaryClimateVariable <: LayerDataset end
struct SecondaryClimateVariable <: LayerDataset end

"""
    BioClim

A list of 19 bioclimatic variables derived from the monthly temperature and
precipitation data. This dataset is provided by `WorldClim` and `CHELSA`, both
of which offer future versions under `CMIP5` and `CMIP6` models.

| Variable | Description                                                |
| ------   | ------                                                     |
| 1        | Annual Mean Temperature                                    |
| 2        | Mean Diurnal Range (Mean of monthly (max temp - min temp)) |
| 3        | Isothermality (BIO2/BIO7) (* 100)                          |
| 4        | Temperature Seasonality (standard deviation *100)          |
| 5        | Max Temperature of Warmest Month                           |
| 6        | Min Temperature of Coldest Month                           |
| 7        | Temperature Annual Range (BIO5-BIO6)                       |
| 8        | Mean Temperature of Wettest Quarter                        |
| 9        | Mean Temperature of Driest Quarter                         |
| 10       | Mean Temperature of Warmest Quarter                        |
| 11       | Mean Temperature of Coldest Quarter                        |
| 12       | Annual Precipitation                                       |
| 13       | Precipitation of Wettest Month                             |
| 14       | Precipitation of Driest Month                              |
| 15       | Precipitation Seasonality (Coefficient of Variation)       |
| 16       | Precipitation of Wettest Quarter                           |
| 17       | Precipitation of Driest Quarter                            |
| 18       | Precipitation of Warmest Quarter                           |
| 19       | Precipitation of Coldest Quarter                           |
"""
struct BioClim <: LayerDataset end

"""
    LandCover

Information on land cover, currently only provided by `EarthEnv`.
"""
struct LandCover <: LayerDataset end

"""
    Elevation

General type for a DEM, currently available through `WorldClim`
"""
struct Elevation <: LayerDataset end

"""
    HabitatHeterogeneity

Information on [habitat heterogeneity](https://www.earthenv.org/texture),
currently only provided by `EarthEnv`.
"""
struct HabitatHeterogeneity <: LayerDataset end

"""
    Topography

Information on habitat topography, currently provided by `EarthEnv`.
"""
struct Topography <: LayerDataset end

# Provider paths
_rasterpath(::Type{TerraClimate}) = "TerraClimate"

# Dataset paths
_rasterpath(::Type{BioClim}) = "BioClim"
_rasterpath(::Type{LandCover}) = "LandCover"
_rasterpath(::Type{HabitatHeterogeneity}) = "HabitatHeterogeneity"
_rasterpath(::Type{Elevation}) = "Elevation"
_rasterpath(::Type{Topography}) = "Topography"
_rasterpath(::Type{PrimaryClimateVariable}) = "Primaries"
_rasterpath(::Type{SecondaryClimateVariable}) = "Secondaries"

