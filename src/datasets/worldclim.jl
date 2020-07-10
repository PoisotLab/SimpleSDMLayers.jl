"""
    worldclim(layers::Vector{T}; resolution::AbstractString="10", path::AbstractString="assets") where {T <: Integer}

Download and prepare WorldClim 2.1 bioclimatic variables, and returns them as an
array of `SimpleSDMPredictor`s. Layers are called by their number, from 1 to 19.
The list of available layers is given in a table below.

The two keywords are `resolution`, which must be a string, and either `2.5`,
`5`, or `10`; and `path`, which refers to the path where the function will look
for the zip and geotiff files.

Internally, this function will download the main zip file for the required
resolution from the WordlClim website, extract it, and parse the required
layers.

It is recommended to *keep* the content of the `path` folder, as it will
eliminate the need to download and/or extract the tiff files. For example,
calling `wordlclim(1:19)` will download and extract everything, and future calls
will be much faster.

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

Original data: https://www.worldclim.org/data/worldclim21.html
"""
function worldclim(layers::Vector{T}; resolution::AbstractString="10", path::AbstractString="assets") where {T <: Integer}
    all(1 .≤ layers .≤ 19) || throw(ArgumentError("The number of the layers must all be between 1 and 19"))
    isdir(path) || mkdir(path)
    resolution ∈ ["2.5", "5", "10"] || throw(ArgumentError("The resolution argument ($(resolution) must be 2.5, 5, or 10"))

    paths = [joinpath(path, "wc2.1_$(resolution)m_bio_$(code).tif") for code in layers]

    missing_files = !all(isfile.(paths))
    zip_file = "bioclim_2.1_$(resolution)m.zip"
    missing_zip_file = !isfile(joinpath(path, zip_file))
    if missing_files
        if missing_zip_file
            @info "Downloading $(zip_file)"
            url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_$(resolution)m_bio.zip"
            r = HTTP.request("GET", url)
            open(joinpath(path, zip_file), "w") do f
                write(f, String(r.body))
            end
        end
        zf = ZipFile.Reader(joinpath(path, zip_file))
        for rf in zf.files
            if joinpath(path, rf.name) in paths
                if !isfile(joinpath(path, rf.name))
                    @info "Reading layer $(rf.name) from archive"
                    write(joinpath(path, rf.name), read(rf))
                end
            end
        end
        close(zf)
    end

    data_layers = geotiff.(paths)

    return SimpleSDMPredictor.(data_layers)

end

"""
    worldclim(layer::T; x...) where {T <: Integer}

Return a single layer from WorldClim 2.0.
"""
worldclim(layer::T; x...) where {T <: Integer} = first(worldclim([layer]; x...))

"""
    worldclim(layers::UnitRange{T}; x...) where {T <: Integer}

Return a range of layers from WorldClim 2.0.
"""
worldclim(layers::UnitRange{T}; x...) where {T <: Integer} = worldclim(collect(layers); x...)
