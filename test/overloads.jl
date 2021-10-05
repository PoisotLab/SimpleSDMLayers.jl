module SSLTestOverloads
using SimpleSDMLayers
using Test

M = rand(Bool, (5, 10))
S = SimpleSDMPredictor(M, 0.0, 1.0, 0.0, 1.0)

for i in eachindex(M)
   @test S[i] == M[i]
end
for i in eachindex(S)
   @test S[i] == M[i]
end

@test typeof(S[0.2, 0.6]) == eltype(M)
@test S[0.2, 0.6] == S[Point(0.2, 0.6)]
@test isnothing(S[1.2, 0.3])
@test isnothing(S[1.2, 1.3])
@test isnothing(S[0.2, 1.3])

C = (left=0.2, bottom=0.5)
@test typeof(clip(S; C...)) == typeof(S)

Y = SimpleSDMResponse(zeros(Float64, (5,5)), 0.0, 1.0, 0.0, 1.0)
Y[0.1,0.1] = 0.2
@test Y[Point(0.1,0.1)] == 0.2

Z = convert(SimpleSDMPredictor, Y)
Y[0.1,0.1] = 4.0
@test Z[0.1,0.1] != Y[0.1,0.1]
@test Z[0.1,0.1] == 0.2

Z = convert(SimpleSDMPredictor, Y)
V = collect(Z)
@test typeof(V) == Vector{Float64}

# typed similar
c2 = similar(Y, Bool)
@test eltype(c2) == Bool

@test eltype(convert(Int64, c2)) == Int64
@test eltype(convert(Float32, c2)) == Float32

# replacement
s1 = SimpleSDMResponse(collect(reshape(1:9, 3, 3)))
@test length(s1) == 9
replace!(s1, 1 => 2, 3 => 2, 9 => nothing)
@test length(s1) == 8
@test s1.grid[1,1] == 2
@test s1.grid[3,1] == 2
@test s1.grid[1,3] == 7
@test isnothing(s1.grid[3,3])

s1 = SimpleSDMPredictor(collect(reshape(1:9, 3, 3)))
@test length(s1) == 9
s2 = replace(s1, 1 => 2, 3 => 2, 9 => nothing)
@test length(s2) == 8
@test s2.grid[1,1] == 2
@test s2.grid[3,1] == 2
@test s2.grid[1,3] == 7

# ==, isequal, hash
l1, l2 = SimpleSDMPredictor(WorldClim, BioClim, 1:2; left = 0.0, right = 10.0, bottom = 0.0, top = 10.0)
l3 = copy(l1)
l4 = similar(l1)
replace!(l4, nothing => NaN)
l5 = SimpleSDMPredictor(replace(l1.grid, nothing => missing), l1)

@test l1 == l1
@test l1 === l1
@test l2 != l1
@test l3 == l1
@test l3 !== l1

@test l4 != l1
@test l4 != l4
@test !isequal(l4, l1)
@test isequal(l4, l4)

@test ismissing(l5 == l1)
@test ismissing(l5 == l5)
@test !isequal(l5, l1)
@test isequal(l5, l5)

@test hash(l1) == hash(l1)
@test hash(l2) != hash(l1)
@test hash(l3) == hash(l1)
@test hash(l4) != hash(l1)
@test hash(l4) == hash(l4)
@test hash(l5) != hash(l1)
@test hash(l5) == hash(l5)

# getindex(layer1, layer2)
l1, l2 = SimpleSDMPredictor(WorldClim, BioClim, 1:2; left = 0.0, right = 10.0, bottom = 0.0, top = 10.0)
l3 = SimpleSDMPredictor(WorldClim, BioClim, 1; left = 5.0, right = 10.0, bottom = 5.0, top = 10.0)
@test stride(l1) == stride(l3)

l4 = clip(l1, l2)
@test l4 == l1

l5 = clip(l1, l3)
@test l5 == l3
@test_throws ArgumentError clip(l3, l1)

end
