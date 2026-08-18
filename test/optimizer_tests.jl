using Test
using RedCarD

@testset "_mutirr" begin
    v = RedCarD._mutirr(5)
    @test length(v) == 5
    @test all(x -> 0.0 <= x < 1.0, v)
end

@testset "errorfind!" begin
    ham = RedCarD.hamiltonian("TFIM", 2, [1, 1], UInt)
    sub = RedCarD.PauliList(["Z-", "-Z"])
    hamcopy = RedCarD.PauliSentence{keytype(ham),ComplexF64}(ham)
    err = RedCarD.errorfind!(hamcopy, sub)
    @test 0.0 <= err <= 1.0
end

@testset "analytic rotosolve cost" begin
    # `_costcoeffs` replaces three explicit `cost` evaluations with the coefficients of
    # c₀ + c₁cos(2θ) + c₂sin(2θ); both must describe the same cost curve, and the minimizer
    # must agree with the three-point formula it replaces.
    Q = 4
    ham = RedCarD.hamiltonian("MFIM", Q, [1, 0.9, 0.4], UInt32)
    elem = RedCarD.PauliSentence(
        RedCarD.PauliList(["Z---", "-Z--", "--ZZ", "XY--"]),
        [0.31, -0.72, 0.5, 1.4],
    )
    elem = RedCarD.PauliSentence{UInt32,ComplexF64}(elem)
    for gs in ["XX--", "-YZ-", "--XY", "ZZZZ"]
        g = RedCarD.UPauli{UInt32}(RedCarD.UPauli(gs))
        c1, c2 = RedCarD._costcoeffs(elem, g, ham)
        c0 = RedCarD.cost(elem, g, 0.0, ham) - c1
        for θ in (0.0, 0.21, pi / 4, pi / 2, -0.83, 1.7)
            @test RedCarD.cost(elem, g, θ, ham) ≈ c0 + c1 * cos(2θ) + c2 * sin(2θ) atol =
                1e-12
        end
        points = [RedCarD.cost(elem, g, i * pi / 4, ham) for i in 0:2]
        @test RedCarD._minanglefind(points) ≈ RedCarD._minanglefind(c1, c2)
    end
end

@testset "optimizer roto minimal" begin
    ham = RedCarD.hamiltonian("TFIM", 4, [1, 1], UInt32)
    d = RedCarD.involutionlessdecomp(RedCarD.PauliList(collect(keys(ham)), 4))
    k = d[:k]
    h = d[:h]

    if (k === nothing) || (length(k) == 0)
        @test true # skip if no k generators
    else
        res = redirect_stdout(devnull) do
            RedCarD.optimizer(ham, h, k; maxiter=0, convergence_tol=1e-6, itertrack=true)
        end
        @test isa(res, NamedTuple)
        @test haskey(res, :H)
        @test haskey(res, :angles)
        @test haskey(res, :error)
        @test haskey(res, :iterations)
        @test res[:iterations] >= 1
        # conjugate back the transformed Hamiltonian and compare to original
        # reverse the generator order and negate the angles when conjugating back
        conjback = RedCarD.ad(res[:H], k, res[:angles])
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
