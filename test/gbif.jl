module SSLTestGBIF
using SimpleSDMLayers
using GBIF
using Test

temperature = worldclim(1)

kingfisher = GBIF.taxon("Megaceryle alcyon", strict=true)

o = GBIF.occurrences(kingfisher)
filter!(GBIF.have_ok_coordinates, o)

# Extract from a single record
for oc in o
    @test typeof(temperature[oc]) <: Number
end

# Modify a value
numboc = SimpleSDMResponse(zeros(Int64, (20, 20)))
for oc in o
    numboc[oc] += 1
end

@info typeof(convert(Matrix, numboc))

@test typeof(convert(Matrix, numboc)) <: Matrix{Integer}

end
