abstract type SimpleSDMSource end

latitudes(::T) where {T <: SimpleSDMSource} = [-90.0, 90.0]
longitudes(::T) where {T <: SimpleSDMSource} = [-180.0, 180.0]

struct WorldClim <: SimpleSDMSource
    resolution::AbstractFloat
    function WorldClim(resolution::AbstractFloat)
        resolution ∈ [2.5, 5.0, 10.0] || throw(ArgumentError("The resolution argument ($(resolution)) must be 2.5, 5, or 10"))
        return new(resolution)
    end
end

WorldClim() = WorldClim(10.0)

struct BioClim <: SimpleSDMSource end

struct EarthEnv <: SimpleSDMSource
    full::Bool
end

EarthEnv() = EarthEnv(false)
latitudes(::EarthEnv) = [-90.0, 90.0]
longitudes(::EarthEnv) = [-156.0, 180.0]
