using Cartan
using Test

@testset "Cartan.jl" begin
    @test product([1, 1, 1, 1], [1, 2, 3, 4]) == ([1, 2, 3, 4], 1, true)
    @test product([2, 2, 2, 2], [1, 2, 3, 4]) == ([2, 1, 4, 3], 1, true)
    @test product([3, 3, 3, 3], [1, 2, 3, 4]) == ([3, 4, 1, 2], 1, true)
    @test product([4, 4, 4, 4], [1, 2, 3, 4]) == ([4, 3, 2, 1], 1, true)
    @test product([1, 2], [3, 4]) == ([3, 3], -1im, false)
    @test product([2, 3, 2], [2, 4, 3]) == ([1, 2, 4], -1, true)
end
