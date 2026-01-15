using Test
using RedCarD

@testset "Pauli string generators" begin
    @test collect(RedCarD.generatex(4, UInt)) == UInt[1, 2, 4, 8]
    @test collect(RedCarD.generatez(4, UInt)) == UInt[16, 32, 64, 128]
    @test collect(RedCarD.generatey(4, UInt)) == UInt[17, 34, 68, 136]
    @test collect(RedCarD.generatexx(4, UInt)) == UInt[3, 6, 12]
    @test collect(RedCarD.generatexx(4, UInt; pbc=true)) == UInt[3, 6, 12, 9]
    @test collect(RedCarD.generatezz(4, UInt)) == UInt[48, 96, 192]
    @test collect(RedCarD.generatezz(4, UInt; pbc=true)) == UInt[48, 96, 192, 144]
    @test collect(RedCarD.generatexy(4, UInt)) == UInt[19, 38, 76]
    @test collect(RedCarD.generateyx(4, UInt)) == UInt[35, 70, 140]
    @test collect(RedCarD.generateyy(4, UInt)) == UInt[51, 102, 204]
end

@testset "Hamiltonians" begin
    ham = hamiltonian("ISING", 4, [2.0], UInt)
    vals = collect(values(ham))
    @test length(vals) == 3
    @test all(isapprox.(abs.(vals), 2.0))
    # wrong couplings
    @test_throws ArgumentError hamiltonian("ISING", 4, [1.0, 2.0], UInt)

    ham = hamiltonian("TFIM", 4, [2.0, 0.5], UInt)
    vals = collect(values(ham))
    @test length(vals) == 7
    @test count(v -> isapprox(abs(v), 2.0), vals) == 3
    @test count(v -> isapprox(abs(v), 1.0), vals) == 4

    ham = hamiltonian("HEISENBERG", 4, [1.5], UInt)
    vals = collect(values(ham))
    @test length(vals) == 9
    @test all(isapprox.(abs.(vals), 1.5))

    # XXZ
    ham = hamiltonian("XXZ", 4, [2.0, -0.8], UInt)
    vals = collect(values(ham))
    @test length(vals) == 9
    @test count(v -> isapprox(abs(v), 2.0), vals) == 6
    @test count(v -> isapprox(abs(v), 1.6), vals) == 3

    # XXZ_EXT: check external field coefficients equal g
    ham = hamiltonian("XXZ_EXT", 4, [2.0, -0.5, 0.25], UInt)
    vals = collect(values(ham))
    @test length(vals) == 13
    @test count(v -> isapprox(abs(v), 0.25), vals) == 4

    ham = hamiltonian("XY", 4, [1.5], UInt)
    vals = collect(values(ham))
    @test length(vals) == 6
    @test all(isapprox.(abs.(vals), 1.5))

    ham = hamiltonian("TFXY", 4, [2.0, 0.5], UInt)
    vals = collect(values(ham))
    @test length(vals) == 10
    @test count(v -> isapprox(abs(v), 2.0), vals) == 6
    @test count(v -> isapprox(abs(v), 1.0), vals) == 4

    ham = hamiltonian("MFIM", 4, [2.0, 0.6, 0.25], UInt)
    vals = collect(values(ham))
    @test length(vals) == 11
    @test count(v -> isapprox(abs(v), 2.0), vals) == 3
    @test count(v -> isapprox(abs(v), 2.0 * 0.6), vals) == 4
    @test count(v -> isapprox(abs(v), 2.0 * 0.25), vals) == 4

    ham = hamiltonian("HIM", 4, [2.0, 0.6, 0.25, 1.0], UInt)
    vals = collect(values(ham))
    @test length(vals) == 11
    @test length(unique(round.(abs.(vals), digits=6))) > 1

    ham = hamiltonian("GN", 4, [0.2, 1.0, 0.7], UInt; nf=1)
    vals = collect(values(ham))
    @test length(vals) == 10
    @test count(v -> isapprox(abs(v), abs(1 - 2 * 0.2)), vals) == 6
    @test count(v -> isapprox(abs(v), 2 * 0.7), vals) == 4

    @test_throws Exception hamiltonian("SYK", 4, [1.0], UInt)

    # unknown model
    @test_throws ArgumentError hamiltonian("UNKNOWN", 4, [1.0], UInt)
end

@testset "PBC models" begin
    ham = hamiltonian("ISING",4,[2.0],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 4
    @test all(isapprox.(abs.(vals), 2.0))

    ham = hamiltonian("XY",4,[1.5],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 8
    @test all(isapprox.(abs.(vals), 1.5))

    ham = hamiltonian("TFIM",4,[2.0,0.5],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 8
    @test count(v->isapprox(abs(v),2.0), vals) == 4
    @test count(v->isapprox(abs(v),1.0), vals) == 4

    ham = hamiltonian("TFXY",4,[2.0,0.5],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 12
    @test count(v->isapprox(abs(v),2.0), vals) == 8
    @test count(v->isapprox(abs(v),1.0), vals) == 4

    ham = hamiltonian("HEISENBERG",4,[1.5],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 12
    @test all(isapprox.(abs.(vals), 1.5))

    ham = hamiltonian("XXZ",4,[2.0,-0.8],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 12
    @test count(v->isapprox(abs(v), 2.0), vals) == 8
    @test count(v->isapprox(abs(v), 1.6), vals) == 4

    ham = hamiltonian("XXZ_EXT",4,[2.0, -0.5, 0.25],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 16
    @test count(v->isapprox(abs(v), 0.25), vals) == 4

    ham = hamiltonian("MFIM",4,[2.0,0.6,0.25],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 12

    ham = hamiltonian("HIM",4,[2.0,0.6,0.25,1.0],UInt; pbc=true)
    vals = collect(values(ham))
    @test length(vals) == 12
end
