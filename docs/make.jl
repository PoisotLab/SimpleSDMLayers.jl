push!(LOAD_PATH, joinpath("..", "src"))

using Documenter, SimpleSDMLayers

makedocs(
    sitename = "Simple SDM Layers"
)

deploydocs(
    repo = "gtihub.com/EcoJulia/SimpleSDMLayers.jl.git"
)
