function iterativeoptimizer(ham::AbstractDict{<:AbstractVector{Int8},<:Real},
    abstrings::AbstractMatrix{Int8},
    symgenerators::AbstractVector{<:AbstractMatrix{Int8}},
    initangles::AbstractVector{<:AbstractVector{<:Real}}=[pi * rand(size(symgen, 2)) for symgen in symgenerators];
    method::AbstractString="roto",
    maxiter::Integer=0,
    convergence_tol::Real=1e-6,
    coeff_tol::Real=0,
    toltype::AbstractString="relerror",
    itertrack::Bool=false,
    timetrack::Bool=false)

    angles = copy(initangles)
    relerror = 1.0
    iter = 0
    stepham = ham
    t = 0
    for i in axes(abstrings, 2)
        ittol = max(1e-4, convergence_tol * 10.0^(1 - i))
        println("Begin optimization for abelian element $i")
        opt = optimizer(stepham, abstrings[:, i:i], symgenerators[i], angles[i], method=method, maxiter=maxiter, convergence_tol=ittol, coeff_tol=coeff_tol, toltype=toltype, itertrack=itertrack, timetrack=timetrack)
        println()
        stepham = opt["H"]
        angles[i] = opt["angles"]
        itertrack && (iter += opt["iterations"])
        timetrack && (t += opt["time"])
    end
    finalham = conjugate(ham, reverse(hcat(symgenerators...), dims=2), -reverse(vcat(angles...)), atol=coeff_tol)
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