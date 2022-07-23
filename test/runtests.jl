using SimpleSDMLayers
using Test

global anyerrors = false

tests = [
    "construction" => "core/construction.jl",
    "lat./lon." => "core/latlon.jl",
    "lat./lon. conversion" => "core/coordconvert.jl",
    "iteration" => "core/iteration.jl",
    "tiling" => "core/tiling.jl",
    "setindex" => "core/setindex.jl",
    "ascii" => "data/ascii.jl",
    "import" => "data/dataread.jl",
    "worldclim" => "data/worldclim.jl",
    "earthenv" => "data/earthenv.jl",
    "chelsa" => "data/chelsa.jl",
    "overloads" => "overloads.jl",
    "clipping" => "operations/clip.jl",
    "rescale" => "operations/rescale.jl",
    "mosaic" => "operations/mosaic.jl",
    "coarsen" => "operations/coarsen.jl",
    "subsetting" => "operations/subsetting.jl",
    "generated" => "generated.jl",
    "plotting" => "extensions/plots.jl",
    "GBIF" => "extensions/gbif.jl",
]

for test in tests
    try
        include(test.second)
        println("\033[1m\033[32m✓\033[0m\t$(test.first)")
    catch e
        global anyerrors = true
        println("\033[1m\033[31m×\033[0m\t$(test.first)")
        println("\033[1m\033[38m→\033[0m\ttest/$(test.second)")
        showerror(stdout, e, backtrace())
        println()
        break
    end
end

if anyerrors
    throw("Tests failed")
end
