module SSLTestCHELSA
using SimpleSDMLayers
using Test

lc1 = SimpleSDMPredictor(CHELSA, BioClim, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(lc1) <: SimpleSDMPredictor

lc1 = SimpleSDMPredictor(CHELSA, BioClim, CanESM2, RCP26, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(lc1) <: SimpleSDMPredictor

lc1 = SimpleSDMPredictor(CHELSA, BioClim, MIROC6, SSP126, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(lc1) <: SimpleSDMPredictor

end
