push!(LOAD_PATH, joinpath("..", "src"))

using Documenter, SimpleSDMLayers
using GBIF
using Statistics

# Set a RasterDataSources path for testing
rdsp = get(ENV, "RASTERDATASOURCES_PATH", "rasterdata")
(rdsp == "rasterdata") && (ENV["RASTERDATASOURCES_PATH"] = rdsp)
isdir(rdsp) || mkdir(rdsp)


makedocs(
    sitename = "Simple SDM Layers",
    modules = [SimpleSDMLayers],
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Types" => "man/types.md",
            "Overloads" => "man/overloads.md",
            "Other operations" => "man/operations.md",
        ],
        "Examples" => [
            "Temperature data" => "examples/temperature.md",
            "DataFrames integration" => "examples/dataframes.md",
            "Importing raster data" => "examples/import.md",
            "Sliding window analysis" => "examples/slidingwindow.md",
            "Landcover data" => "examples/landcover.md",
            "Landcover consensus" => "examples/consensus.md"
        ],
        "Building SDMs" => [
            "GBIF integration" => "sdm/gbif.md",
            "BIOCLIM from scratch" => "sdm/bioclim.md"
        ]
    ],
    checkdocs=:all,
    strict=true
)

run(`find . -type f -size +40M -delete`)

deploydocs(
    repo = "github.com/EcoJulia/SimpleSDMLayers.jl.git",
    push_preview = true
)
