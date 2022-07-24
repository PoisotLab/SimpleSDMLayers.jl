module SSLTestCHELSA
using SimpleSDMLayers
using Test

lc1 = SimpleSDMPredictor(CHELSA, BioClim, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
@test typeof(lc1) <: SimpleSDMPredictor

for mod in [GFDLESM4, IPSLCM6ALR, MPIESM12HR, MRIESM20, UKESM10LL]
    lc1 = SimpleSDMPredictor(CHELSA, BioClim, mod, SSP126, 1; left=-5.0, right=7.0, bottom=30.0, top=45.0)
    @test typeof(lc1) <: SimpleSDMPredictor
end

end
