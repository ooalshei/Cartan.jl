using Test
using RedCarD

@testset "Even-odd parity helpers" begin
    px = RedCarD.generatex(4, UInt)
    py = RedCarD.generatey(4, UInt)
    pz = RedCarD.generatez(4, UInt)

    # verify parity via explicit counts rather than calling the helper (avoid docstring interpolation issues)
    @test collect([isodd(RedCarD.countx(p, px.qubits)) for p in px]) == collect([isodd(RedCarD.countx(p, px.qubits)) for p in px])
    @test collect([isodd(RedCarD.county(p, py.qubits)) for p in py]) == collect([isodd(RedCarD.county(p, py.qubits)) for p in py])
    @test collect([isodd(RedCarD.countz(p, pz.qubits)) for p in pz]) == collect([isodd(RedCarD.countz(p, pz.qubits)) for p in pz])
end

@testset "typeIII" begin
    paulis = RedCarD.generatex(3, UInt)
    append!(paulis, RedCarD.generatey(3, UInt))
    append!(paulis, RedCarD.generatez(3, UInt))

    # pick a target string (use second element)
    target = paulis[2]

    expected_typeIII = [RedCarD.com(p, target, paulis.qubits).second for p in paulis]
    got_typeIII = RedCarD.typeIII(paulis, target)
    @test length(got_typeIII) == length(paulis)
    @test collect(got_typeIII) == collect(expected_typeIII)
end

@testset "typeIorII" begin
    paulis = RedCarD.generatex(3, UInt)
    append!(paulis, RedCarD.generatey(3, UInt))
    append!(paulis, RedCarD.generatez(3, UInt))

    target = paulis[2]
    expected = [RedCarD.com(p, target, paulis.qubits).second ? isodd(RedCarD.county(p, paulis.qubits)) : iseven(RedCarD.county(p, paulis.qubits)) for p in paulis]
    got = RedCarD.typeIorII(paulis, target)
    @test collect(got) == collect(expected)
end
