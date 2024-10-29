function iterativeoptimizer(ham::Dict{Vector{Int8},Float64},
    abstrings::Matrix{Int8},
    symgenerators::Vector{Matrix{Int8}},
    initangles::Vector{Vector{Float64}}=[pi * rand(size(symgen, 2)) for symgen in symgenerators];
    method::String="roto",
    maxiter::Int=0,
    mintol::Float64=1e-6,
    tol::Float64=0.0,
    toltype::String="relerror")::Dict{String,Union{Dict{Vector{Int8},Float64},Vector{Vector{Float64}}}}

    angles = copy(initangles)
    stepham = ham
    for i in axes(abstrings, 2)
        println("Begin optimization for abelian element $i")
        opt = optimizer(stepham, abstrings[:, i:i], symgenerators[i], initangles[i], method=method, maxiter=maxiter, mintol=mintol, tol=tol, toltype=toltype)
        println()
        stepham = opt["H"]
        angles[i] = opt["angles"]
    end
    finalham = conjugate(ham, reverse(hcat(symgenerators...), dims=2), -reverse(vcat(angles...)), tol=tol)
    relerror = errorfind!(finalham, abstrings)

    println("Optimization complete. Combined relative error: $(sqrt(relerror))")
    return Dict("H" => finalham, "angles" => angles)
end
