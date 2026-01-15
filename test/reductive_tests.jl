using Test
using RedCarD

@testset "subalgred" begin
    h = RedCarD.PauliList(["Z--", "-Z-", "ZZ-", "--Z", "Z-Z", "-ZZ", "ZZZ"])
    b = RedCarD.subalgred(h)
    expected = RedCarD.PauliList(["Z--", "-Z-", "--Z"])

    @test length(b) == length(expected)
    for s in expected
        @test s in b
    end
end

@testset "fragmentedsubspaces" begin
    # build a small Hamiltonian and use involutionlessdecomp to obtain k and h
    ham = RedCarD.hamiltonian("TFIM", 4, [1, 1], UInt32)
    d = RedCarD.involutionlessdecomp(RedCarD.PauliList(collect(keys(ham)), 4))
    k = d[:k]
    h = d[:h]
    parts = RedCarD.fragmentedsubspaces(k, h)

    @test isa(parts, Vector)
    @test length(parts) == length(h)
    for p in parts
        @test isa(p, RedCarD.PauliList)
    end
end

@testset "cleangenerators!" begin
    # prepare symgenerators with empty PauliLists interleaved
    empty1 = RedCarD.PauliList{UInt,1}(undef, 0)
    nonempty = RedCarD.PauliList(["Z-"])
    empty2 = RedCarD.PauliList{UInt,1}(undef, 0)
    symgens = [empty1, nonempty, empty2]
    bstrings = RedCarD.PauliList(["Z-", "X-", "Y-"])

    RedCarD.cleangenerators!(symgens, bstrings)
    @test length(symgens) == 1
    @test length(bstrings) == 1
    @test isa(symgens[1], RedCarD.PauliList)
end
