function SimpleSDMLayers.resample(x) 
    GeoData.resample(x)
end

"""
    Code from GeoData.jl

    warp(A::AbstractGeoArray, flags::Dict)

function warp(A::AbstractGeoArray, flags::Dict)
    odims = otherdims(A, (X, Y, Band))
    if length(odims) > 0
        # Handle dimensions other than X, Y, Band
        slices = slice(A, odims)
        warped = map(A -> _warp(A, flags), slices)
        return combine(warped, odims)
    else
        return _warp(A, flags)
    end
end
warp(st::AbstractGeoStack, flags::Dict) = map(A -> warp(A, flags), st)

function _warp(A::AbstractGeoArray, flags::Dict)
    flagvect = reduce([flags...]; init=[]) do acc, (key, val)
        append!(acc, String[_asflag(key), _stringvect(val)...])
    end
    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], flagvect) do warped
            _maybe_permute_from_gdal(read(GeoArray(warped)), dims(A))
        end
    end
end

"""