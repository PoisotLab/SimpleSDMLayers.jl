module SSLTestWorldclim
using SimpleSDMLayers
using Test

temp = SimpleSDMPredictor(WorldClim, BioClim, 1)
@test size(temp) == (1080, 2160)
@test round(first(stride(temp)); digits=2) ≈ round(last(stride(temp)); digits=2)
@test length(longitudes(temp)) == 2160
@test length(latitudes(temp)) == 1080

wc3 = SimpleSDMPredictor(WorldClim, BioClim, 3)
@test typeof(wc3) <: SimpleSDMPredictor

future = SimpleSDMPredictor(WorldClim, BioClim, MIROC6, SSP126, 1)
@test typeof(future) <: SimpleSDMPredictor

future = SimpleSDMPredictor(WorldClim, BioClim, CanESM5, SSP370, 2; year="2061-2080")
@test typeof(future) <: SimpleSDMPredictor

end
