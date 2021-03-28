module SSLTestDataRead
using SimpleSDMLayers
using Test

l = SimpleSDMPredictor(WorldClim, BioClim, 1)
f = tempname()
geotiff(f, l)
mp = geotiff(SimpleSDMResponse, f)

@test typeof(mp) <: SimpleSDMResponse
@test size(mp) == size(l)

end
