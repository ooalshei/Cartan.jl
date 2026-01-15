using Test
using RedCarD

@testset "reductive_optimizer basic" begin
    ham = RedCarD.hamiltonian("TFIM", 4, [1, 1], UInt32)
    paulis = RedCarD.PauliList(collect(keys(ham)), 4)
    d = RedCarD.involutionlessdecomp(paulis)
    bstrings = d[:h]
    k = d[:k]

    if (k === nothing) || (length(k) == 0) || (length(bstrings) == 0)
        @test true # skip when there is nothing to optimize
    else
        symgens = RedCarD.fragmentedsubspaces(k, bstrings)
        @test isa(symgens, Vector)
        @test all(x -> isa(x, RedCarD.PauliList), symgens)

        # run a short reductive optimization with tracking
        res = redirect_stdout(devnull) do
            RedCarD.reductive_optimizer(
                ham,
                bstrings,
                symgens;
                maxiter=0,
                convergence_tol=1e-6,
                itertrack=true,
            )
        end
        @test isa(res, Dict)
        @test haskey(res, :H)
        @test haskey(res, :angles)
        @test haskey(res, :error)
        @test haskey(res, :iterations)
        @test 0.0 <= res[:error] <= 1.0
        @test res[:iterations] >= 0
        # conjugate back the transformed Hamiltonian using concatenated generators and angles
        gens_all = vcat(symgens...)
        angles_all = vcat(res[:angles]...)
        conjback = RedCarD.ad(res[:H], gens_all, angles_all)
        diff = 0.0
        full = 0.0
        for k in union(collect(keys(ham)), collect(keys(conjback)))
            v0 = get(ham, k, 0.0 + 0.0im)
            v1 = get(conjback, k, 0.0 + 0.0im)
            diff += abs2(v0 - v1)
            full += abs2(v0)
        end
        relerr = sqrt(diff / max(full, eps()))
        @test relerr <= 1e-6
    end
end
