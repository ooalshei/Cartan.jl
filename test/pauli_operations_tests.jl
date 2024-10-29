using Cartan
using Test

sentenceA = Dict(Int8[1, 1, 2] => 1.0, Int8[1, 3, 4] => -0.7, Int8[4, 4, 2] => 0.3)
    sentenceB = Dict(Int8[1, 1, 2] => 3.5, Int8[1, 2, 3] => 1.2, Int8[4, 4, 2] => -0.295)
    sentenceC = Dict(Int8[1, 1, 2] => 3.5, Int8[1, 2, 3] => 1.2, Int8[4, 4, 2] => 0.295)
    sentenceAplusBplusC = Dict(Int8[1, 1, 2] => 8.0, Int8[1, 3, 4] => -0.7, Int8[1, 2, 3] => 2.4, Int8[4, 4, 2] => 0.3)
    sentenceCminusBminusA = Dict(Int8[1, 1, 2] => -1.0, Int8[1, 3, 4] => 0.7, Int8[4, 4, 2] => 0.29)
    sentenceAB = Dict(Int8[1, 1, 1] => 3.4115, Int8[1, 2, 4] => 1.2im, Int8[4, 4, 1] => 0.755, Int8[1, 3, 3] => -2.45im, Int8[1, 4, 2] => 0.84, Int8[4, 2, 3] => -0.2605, Int8[4, 3, 4] => -0.36)

@testset "Cartan.jl" begin
    @test pauliprod(Int8[1, 1, 1, 1], Int8[1, 2, 3, 4]) == ([1, 2, 3, 4], 1, true)
    @test pauliprod(Int8[2, 2, 2, 2], Int8[1, 2, 3, 4]) == ([2, 1, 4, 3], 1, true)
    @test pauliprod(Int8[3, 3, 3, 3], Int8[1, 2, 3, 4]) == ([3, 4, 1, 2], 1, true)
    @test pauliprod(Int8[4, 4, 4, 4], Int8[1, 2, 3, 4]) == ([4, 3, 2, 1], 1, true)
    @test pauliprod(Int8[1, 2], Int8[3, 4]) == ([3, 3], -1im, false)
    @test pauliprod(Int8[2, 3, 2], Int8[2, 4, 3]) == ([1, 2, 4], -1, true)
    @test paulisum(sentenceA, sentenceB, sentenceC) == sentenceAplusBplusC
    @test paulidiff(sentenceC, sentenceB, sentenceA) == sentenceCminusBminusA
    
    # prodAB = pauliprod(sentenceA, sentenceB; tol=0.001)
    # testarray = copy(sentenceAB)
    # for (key, value) in prodAB
    #     abs(value - testarray[key]) <= 0.001 && pop!(testarray, key)
    # end
    # @test length(testarray) == 0
end
@testset "Cartan.jl" begin
    ABCconjA = Dict(Int8[1, 1, 2] => 1.3597371432019283, Int8[1, 3, 3] => -3.917242500465444, Int8[1, 3, 4] => -6.877232444373558, Int8[1, 2, 3] => -0.9987524077131418, Int8[1, 2, 4] => -2.182313824381636, Int8[4, 4, 2] => 0.3)
    angles = [1.0, -0.7, 0.3]
    generators = Int8[[1, 1, 2] [1, 3, 4] [4, 4, 2]]
    @test conjugate(sentenceAplusBplusC, generators, angles) == ABCconjA
end
