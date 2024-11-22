function iterativeoptimizer(ham::Dict{Vector{Int8},Float64},
    abstrings::Matrix{Int8},
    symgenerators::Vector{Matrix{Int8}},
    initangles::Vector{Vector{Float64}}=[pi * rand(size(symgen, 2)) for symgen in symgenerators];
    method::String="roto",
    maxiter::Integer=0,
    mintol::Float64=1e-6,
    tol::Float64=0.0,
    toltype::String="relerror",
    itertrack::Bool=false,
    timetrack::Bool=false)::Dict{String,Union{Dict{Vector{Int8},Float64},Vector{Vector{Float64}},Integer,Float64}}

    angles = copy(initangles)
    relerror = 1.0
    iter = 0
    stepham = ham
    t = 0
    for i in axes(abstrings, 2)
        ittol = max(1e-4, mintol * 10.0^(1 - i))
        println("Begin optimization for abelian element $i")
        opt = optimizer(stepham, abstrings[:, i:i], symgenerators[i], angles[i], method=method, maxiter=maxiter, mintol=ittol, tol=tol, toltype=toltype, itertrack=itertrack, timetrack=timetrack)
        println()
        stepham = opt["H"]
        angles[i] = opt["angles"]
        itertrack && (iter += opt["iterations"])
        timetrack && (t += opt["time"])
    end
    finalham = conjugate(ham, reverse(hcat(symgenerators...), dims=2), -reverse(vcat(angles...)), tol=tol)
    relerror = errorfind!(finalham, abstrings)
    println("Combined relative error: $(sqrt(relerror))")

    println("Optimization complete. Combined relative error: $(sqrt(relerror))")
    if itertrack
        timetrack ? (return Dict("H" => finalham, "angles" => angles, "iterations" => iter, "time" => t)) : (return Dict("H" => finalham, "angles" => angles, "iterations" => iter))
    else
        timetrack ? (return Dict("H" => finalham, "angles" => angles, "time" => t)) :
        (return Dict("H" => finalham, "angles" => angles))
    end
end
