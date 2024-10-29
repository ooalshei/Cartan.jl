using Cartan
using SafeTestsets

@testset "Cartan.jl" begin
    ham = hamiltonian("TFIM", 4, [1.0, 1.0])
    cd = involutionlessdecomp(ham[1])
    hamdict = Dict(ham[1][:, i] => ham[2][i] for i in eachindex(ham[2]))
    initangles = zeros(Float64, size(cd["k"], 2))
    H = Dict{Vector{Int8}, Float64}([1, 1, 1, 4] => -1.0000000000000469, [1, 4, 1, 1] => -1.879385241566419, [4, 1, 1, 1] => -0.347296355333871, [1, 1, 4, 1] => -1.5320888862432944)
    angles = [3.4608460344327567, 2.5233281491570123, 3.5116000402567646, 2.731011961453448, 2.9170218274717, 3.273026982517882, 3.3714923713665694, 2.7332931487733614, 2.8914873116346285, 3.22244825315624, 2.6006886819233324, 3.4307291707529775]

    touch("temp.txt")
    open("temp.txt", "w") do io
        redirect_stdout(io) do
            opt = optimizer(hamdict, cd["h"], cd["k"], initangles, method="roto", maxiter=0, mintol=1e-6, tol=0.0, toltype="relerror")
            @test (paulidiff(opt["H"], H, tol=1e-3) == Dict{Vector{Int8}, Float64}())
        end
    end
    rm("temp.txt")
end

@testset "Cartan.jl" begin
    model = "TFIM"
    sites = 4
    couplings = [1.0, 1.0]
    ham = hamiltonian(model, sites, couplings)
    cd = involutionlessdecomp(ham[1])
    hamdict = Dict(ham[1][:, i] => ham[2][i] for i in eachindex(ham[2]))
    abstrings = subalgred(cd["h"])
    symstrings = symsubspaces(cd["k"], abstrings)
    abstrings = cleangenerators!(symstrings, abstrings)
    initangles = [zeros(Float64, size(symgen, 2)) for symgen in symstrings]
    H = Dict{Vector{Int8}, Float64}([1, 1, 1, 4] => -0.34728712503478937, [1, 4, 1, 1] => -1.5323146339960143, [4, 1, 1, 1] => -1.8791252998862664, [1, 1, 4, 1] => -1.0000175642122628)
    angles = [[3.6378540651442277, 2.544258885895675, 3.54829993220899, 2.8242213350616465, 3.2566110914590922, 2.9201627848578173], [3.6803679303270815, 2.49364614340446, 3.548299933474752, 2.9137272313341613], [3.502246796204928, 2.5442588855361263]]

    touch("temp.txt")
    open("temp.txt", "w") do io
        redirect_stdout(io) do
            opt = iterativeoptimizer(hamdict, abstrings, symstrings, initangles, method="roto", maxiter=0, mintol=1e-10, tol=0.0, toltype="relerror")
            @test (paulidiff(opt["H"], H, tol=1e-3) == Dict{Vector{Int8}, Float64}())
        end
    end
    rm("temp.txt")
end