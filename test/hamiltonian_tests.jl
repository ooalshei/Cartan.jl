using Cartan
using Test

@testset "Cartan.jl" begin
    @test Cartan.generatexx(4, pbc=true) == Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2]]
    @test Cartan.generateyy(4, pbc=true) == Int8[[3, 3, 1, 1] [1, 3, 3, 1] [1, 1, 3, 3] [3, 1, 1, 3]]
    @test Cartan.generatezz(4, pbc=true) == Int8[[4, 4, 1, 1] [1, 4, 4, 1] [1, 1, 4, 4] [4, 1, 1, 4]]
    @test Cartan.generatex(4) == Int8[[2, 1, 1, 1] [1, 2, 1, 1] [1, 1, 2, 1] [1, 1, 1, 2]]
    @test Cartan.generatey(4) == Int8[[3, 1, 1, 1] [1, 3, 1, 1] [1, 1, 3, 1] [1, 1, 1, 3]]
    @test Cartan.generatez(4) == Int8[[4, 1, 1, 1] [1, 4, 1, 1] [1, 1, 4, 1] [1, 1, 1, 4]]
end

@testset "Cartan.jl" begin
    @test hamiltonian("Ising", 4, [0.5], pbc=true) == Tuple([Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2]], [-0.5, -0.5, -0.5, -0.5]])
    @test hamiltonian("XY", 4, [0.5], pbc=true) == Tuple([Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2] [3, 3, 1, 1] [1, 3, 3, 1] [1, 1, 3, 3] [3, 1, 1, 3]], [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5]])
    @test hamiltonian("TFIM", 4, [0.5, 0.0], pbc=true) == Tuple([Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2]], [-0.5, -0.5, -0.5, -0.5]])
    @test hamiltonian("TFIM", 4, [3.0, 0.5], pbc=true) == Tuple([Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2] [4, 1, 1, 1] [1, 4, 1, 1] [1, 1, 4, 1] [1, 1, 1, 4]], [-3.0, -3.0, -3.0, -3.0, -1.5, -1.5, -1.5, -1.5]])
    @test hamiltonian("TFXY", 4, [0.5, 0.0], pbc=true) == Tuple([Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2] [3, 3, 1, 1] [1, 3, 3, 1] [1, 1, 3, 3] [3, 1, 1, 3]], [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5]])
    @test hamiltonian("TFXY", 4, [3.0, 0.5], pbc=true) == Tuple([Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2] [3, 3, 1, 1] [1, 3, 3, 1] [1, 1, 3, 3] [3, 1, 1, 3] [4, 1, 1, 1] [1, 4, 1, 1] [1, 1, 4, 1] [1, 1, 1, 4]], [-3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -1.5, -1.5, -1.5, -1.5]])
    @test hamiltonian("Heisenberg", 4, [2.0], pbc=true) == Tuple([Int8[[2, 2, 1, 1] [1, 2, 2, 1] [1, 1, 2, 2] [2, 1, 1, 2] [3, 3, 1, 1] [1, 3, 3, 1] [1, 1, 3, 3] [3, 1, 1, 3] [4, 4, 1, 1] [1, 4, 4, 1] [1, 1, 4, 4] [4, 1, 1, 4]], [-2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0]]) 
end