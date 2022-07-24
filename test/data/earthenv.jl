module SSLTestLandCover
using SimpleSDMLayers
using Test

lc1 = SimpleSDMPredictor(EarthEnv, LandCover, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(lc1) <: SimpleSDMPredictor

lc1 = SimpleSDMPredictor(EarthEnv, LandCover, 1; full=true, left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(lc1) <: SimpleSDMPredictor

he1 = SimpleSDMPredictor(EarthEnv, HabitatHeterogeneity, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(he1) <: SimpleSDMPredictor

to1 = SimpleSDMPredictor(EarthEnv, Topography, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(to1) <: SimpleSDMPredictor

end
