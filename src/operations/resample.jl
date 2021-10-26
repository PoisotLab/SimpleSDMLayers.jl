#=

The code in this file is adaoted from GeoData.jl
(https://github.com/rafaqz/GeoData.jl), see this issue
(https://github.com/EcoJulia/SimpleSDMLayers.jl/issues/131)

See the following GeoData.jl license (MIT)

----------------------------------------------------------------------
Copyright (c) 2019 Rafael Schouten and Michael Krabbe Borregaard

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
----------------------------------------------------------------------
=#

function resample(A::T) end

"""
Code from GeoData.jl

`method`: A `Symbol` or `String` specifying the method to use for
resampling. Defaults to `:near` (nearest neighbor resampling). See
[resampling
method](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r)
in the gdalwarp docs for a complete list of possible values.

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
