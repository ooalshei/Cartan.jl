using Test
using RedCarD

@testset "Involutionless Cartan decomposition - basic" begin
    paulis = RedCarD.generatex(1, UInt)
    append!(paulis, RedCarD.generatez(1, UInt))

    comp = RedCarD.involutionlessdecomp(paulis)
    @test isa(comp, NamedTuple)

    g = comp[:g]
    h = comp[:h]
    k = comp[:k]
    m = comp[:m]

    @test isa(g, RedCarD.PauliList)
    @test isa(h, RedCarD.PauliList)

    # original generators should be in the generated algebra
    @test all(p -> p in g, paulis)

    # k and m are either PauliLists or `nothing` (no decomposition)
    @test (k === nothing) || isa(k, RedCarD.PauliList)
    @test (m === nothing) || isa(m, RedCarD.PauliList)

    # if k and m are present they should be disjoint and partition g
    if isa(k, RedCarD.PauliList) && isa(m, RedCarD.PauliList)
        @test length(k) + length(m) == length(g)
        for a in k
            @test !(a in m)
        end
        for b in m
            @test !(b in k)
        end
        # h should be subset of m
        for hh in h
            @test hh in m
        end
    else
        # otherwise h should be subset of g
        for hh in h
            @test hh in g
        end
    end
end
