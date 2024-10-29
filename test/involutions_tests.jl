using Cartan
using Test

@testset "Cartan.jl" begin
    strings = Int8[[1, 1, 2] [1, 3, 4] [4, 4, 2]]
    @test evenoddx(strings) == [true, false, true]
    @test evenoddy(strings) == [false, true, false]
    @test evenoddz(strings) == [false, true, false]
end