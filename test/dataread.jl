module SSLTestDataRead
using SimpleSDMLayers
using Test
using RasterDataSources

file = joinpath(dirname(pathof(SimpleSDMLayers)), "..", "data", "connectivity.tiff")
struct MyConnectivityMap <: RasterDataSources.RasterDataSource end
SimpleSDMLayers.latitudes(::Type{MyConnectivityMap}) = (-10.0, 10.0)
SimpleSDMLayers.longitudes(::Type{MyConnectivityMap}) = (-20.0, 20.0)
mp = SimpleSDMLayers.geotiff(SimpleSDMResponse, MyConnectivityMap, file)

@test typeof(mp) <: SimpleSDMResponse
@test size(mp) == (1255, 1206)

end
