function cost(generators::AbstractVector{<:UPauli},
    angles::AbstractVector{Float64},
    subalgelem::PauliSentence,
    ham::PauliSentence;
    atol::Real=0)::Float64
    raw"""
    Calculates the cost function :math:`\mathrm{Tr}(KvK^\dag, \mathcal{H})`.
    """
    sentence = ad(subalgelem, generators, angles, atol=atol)
    result::Float64 = 0.0
    for (key, value) in ham
        if haskey(sentence, key)
            result += (-1)^county(key, ham.qubits) * value * sentence[key]
        end
    end
    return result
end

function _cost!(ret::AbstractVector{<:Real},
    generator::UPauli,
    partialelem::PauliSentence,
    partialham::PauliSentence;
    atol::Real=0)
    raw"""
    Calculates the cost function :math:`\mathrm{Tr}(KvK^\dag, \mathcal{H})` at 0, pi/4, pi/2.
    """

    @sync begin
        Threads.@spawn ret[1] = cost([generator], [0.0], partialelem, partialham, atol=atol)
        Threads.@spawn ret[2] = cost([generator], [pi / 4], partialelem, partialham, atol=atol)
        Threads.@spawn ret[3] = cost([generator], [pi / 2], partialelem, partialham, atol=atol)
    end
    nothing
end

function _rotostep!(angles::AbstractVector{<:Real},
    points::AbstractVector{<:Real},
    partialelem::PauliSentence,
    ham::PauliSentence,
    generators::AbstractVector{<:UPauli};
    atol::Real)

    task = Threads.@spawn ham
    for i in eachindex(angles)
        partialelem = ad(partialelem, generators[i], -angles[i], atol=atol)
        partialham = fetch(task)
        _cost!(points, generators[i], partialelem, partialham, atol=atol)
        angles[i] = _minanglefind(points)
        task = Threads.@spawn ad(partialham, generators[i], -angles[i], atol=atol)
        # cosines[i] = cos(2 * angles[i])
        # sines[i] = -sin(2 * angles[i])
    end
    return fetch(task)
end

function errorfind!(ham::PauliSentence, subalgebra::AbstractVector{<:UPauli{T,Q}})::Float64 where {T,Q}
    """
    errorfind
    ----------
    Finds the error of the Hamiltonian after the transformation.
    """
    errornorm = 0.0
    fullnorm = 0.0
    for (key, value) in ham
        fullnorm += abs2(value)
        unique!(com.(subalgebra, [UPauli{T,Q}(key)])) == [UPauli{T,Q}(0)] || (errornorm += abs2(value); delete!(ham, key))
    end
    return errornorm / fullnorm
end

function optimizer(ham::PauliSentence,
    subalgebra::AbstractVector{<:UPauli},
    generators::AbstractVector{<:UPauli},
    initangles::AbstractVector{<:Real}=pi * rand(length(generators));
    method::AbstractString="roto",
    maxiter::Integer=0,
    convergence_tol::Real=1e-6,
    coeff_tol::Real=0,
    toltype::AbstractString="relerror",
    itertrack::Bool=false,
    timetrack::Bool=false)

    length(initangles) == length(generators) || throw(ArgumentError("Incorrect number of initial angles. Expected $(length(generators)), got $(length(initangles))."))

    irr = _mutirr(length(subalgebra))
    subalgelem = PauliSentence{UInt, Float64}(subalgebra, irr)
    if method == "roto"
        angles = copy(initangles)
        # cosines = cos.(2 .* angles)
        # sines = .-sin.(2 .* angles)
        points = Vector{Float64}(undef, 3)
        errorcache = 1.0

        iter = 0
        t = time()
        while true
            iter += 1
            # println("Begin iteration $iter...")
            partialelem = ad(subalgelem, generators, angles, atol=coeff_tol)
            transformedham = _rotostep!(angles, points, partialelem, ham, generators, atol=coeff_tol)
            if (iter % 10 == 0) | (iter == maxiter)
                if toltype == "relerror"
                    # transformedham = conjugate(ham, reverse(generators, dims=2), -reverse(angles), atol=coeff_tol)
                    relerror = errorfind!(transformedham, subalgebra)
                    if (relerror <= convergence_tol^2) | (iter == maxiter)
                        iter == maxiter ? println("Max iterations reached.") : println("Converged in $iter iterations.")
                        println("Final relative error: $(sqrt(relerror))")

                        if itertrack
                            timetrack ? (return Dict("H" => transformedham, "angles" => angles, "iterations" => iter, "time" => time() - t)) : (return Dict("H" => transformedham, "angles" => angles, "iterations" => iter))
                        else
                            timetrack ? (return Dict("H" => transformedham, "angles" => angles, "time" => time() - t)) : (return Dict("H" => transformedham, "angles" => angles))
                        end
                    elseif (sqrt(errorcache) - sqrt(relerror)) <= 0
                        println("Relative error after $iter iterations: $(sqrt(relerror))")
                        println("Convergence too slow. Starting over.")
                        println()
                        iter = 0
                        angles = pi * rand(length(generators))
                        t = time()
                        errorcache = 1.0
                        # angles +=  convergence_tol * randn(size(angles))

                    else
                        errorcache = relerror
                        iter % 50 == 0 && println("Relative error after $iter iterations: $(sqrt(relerror))")
                        println()
                    end
                end
            end
        end
    end
end