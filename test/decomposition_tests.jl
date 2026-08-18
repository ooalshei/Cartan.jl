using Test
using RedCarD

@testset "DLA (dynamical Lie algebra)" begin
    paulis = RedCarD.generatex(1, UInt)
    append!(paulis, RedCarD.generatez(1, UInt))
    closure = RedCarD.dla(paulis)

    # original generators present
    @test all(p -> p in closure, paulis)

    # commutator of X and Z should generate Y
    y = RedCarD.generatey(1, UInt)[1]
    @test y in closure
end

@testset "Subalgebra finder" begin
    paulis = RedCarD.generatez(3, UInt) # all Z strings commute
    sub = RedCarD.subalgfind(paulis)
    # sub should be abelian: all pairs commute
    for i in eachindex(sub)
        for j in (i+1):length(sub)
            @test RedCarD.com(sub[i], sub[j], sub.qubits).second
        end
    end
end

@testset "Cartan decomposition" begin
    paulis = RedCarD.generatex(2, UInt)
    append!(paulis, RedCarD.generatez(2, UInt))

    comp = RedCarD.cartandecomp(paulis, RedCarD.evenoddx)
    @test isa(comp, NamedTuple)
    k = comp[:k]
    m = comp[:m]
    h = comp[:h]

    # partition: k and m together equal original set and are disjoint
    @test length(k) + length(m) == length(paulis)
    for a in k
        @test !(a in m)
    end
    for b in m
        @test !(b in k)
    end

    # h is subset of m
    for hh in h
        @test hh in m
    end

    # h should be abelian
    for i in eachindex(h)
        for j in (i+1):length(h)
            @test RedCarD.com(h[i], h[j], h.qubits).second
        end
    end
end
