raw"""
cartan_optimizer
----------------
This module finds the transformation that places the Hamiltonian into the Cartan subalgebra. For a Cartan Hamiltonian,
we use the KHK theorem which states that for any given term in :math:`\mathfrak{m}`, we can write
:math:`\mathfrak{m} = K \mathfrak{h} K^\dagger` where :math:`K \in \mathrm{e}^{i\mathfrak{k}}`. If we pick a
:math:`v \in \mathfrak{h}` such that :math:`\mathrm{e}^{iv}` is dense in :math:`\mathrm{e}^{i\mathfrak{h}}`, then
optimizing the function :math:`\langle K v K^\dagger, H\rangle` would find the KHK decomposition.

For a non-Cartan Hamiltonian, we can use the same algorithm but we use :math:`G \in \mathrm{e}^{i\mathfrak{g}}` instead
of :math:`K`.
"""
function _mutirr(len::Integer; seed::Irrational=pi)::Vector{Float64}
    """
    _mutirr
    -------
    Generate a random vector of length `len` of mutually irrational numbers.
    """
    result = Vector{Float64}(undef, len)
    result[1] = seed % 1
    for i in 2:len
        result[i] = (seed * result[i-1]) % 1
    end
    return result
end

# @. cosine(x, p) = p[1] * cos(2 * x + p[2]) + p[3] 

function _minanglefind(points::AbstractVector{<:Real})::Float64
    """
    _minanglefind
    ----------
    Finds the angle that minimizes a sine curve that passes through three points at 0, pi/2, pi.
    """
    # cons = 0.5 * (points[1] + points[3])
    # return (-0.5 * atan(points[1] - cons, points[2] - cons) - pi / 4)
    # a = sqrt((points[1] - cons)^2 + (points[2] - cons)^2)
    # points[1] > cons ? (delta = asin((cons - points[2]) / a)) : (delta = pi - asin((cons - points[2]) / a))
    # return 0.5 * (pi - delta)
    # fit = curve_fit(cosine, [0, pi /4, pi / 2], points, [1.0, 0.0, 0.0])
    # d = coef(fit)[2]
    # coef(fit)[1] < 0 ? (return -d / 2) : (return (pi - d) / 2)
    return 0.5 * atan(2 * points[2] - (points[1] + points[3]), (points[1] - points[3])) + pi / 2

end

function cost(generators::AbstractMatrix{Int8},
    angles::AbstractVector{Float64},
    subalgelem::AbstractDict{<:AbstractVector{Int8},<:Real},
    ham::AbstractDict{<:AbstractVector{Int8},<:Real};
    atol::Real=0)::Float64
    raw"""
    Calculates the cost function :math:`\mathrm{Tr}(KvK^\dag, \mathcal{H})`.
    """
    sentence = conjugate(subalgelem, generators, angles, atol=atol)
    costval = 0.0
    for (key, value) in ham
        key in keys(sentence) && (costval += value * sentence[key])
    end
    return costval
end

function _cost!(ret::AbstractVector{<:Real},
    generator::AbstractMatrix{Int8},
    partialelem::AbstractDict{<:AbstractVector{Int8},<:Real},
    partialham::AbstractDict{<:AbstractVector{Int8},<:Real};
    atol::Real=0)
    raw"""
    Calculates the cost function :math:`\mathrm{Tr}(KvK^\dag, \mathcal{H})` at 0, pi/4, pi/2.
    """

    @sync begin
        Threads.@spawn ret[1] = cost(zeros(Int8, 0, 0), zeros(Float64, 0), partialelem, partialham, atol=atol)
        Threads.@spawn ret[2] = cost(generator, [pi / 4], partialelem, partialham, atol=atol)
        Threads.@spawn ret[3] = cost(generator, [pi / 2], partialelem, partialham, atol=atol)
    end
    nothing
end

function _rotostep!(angles::AbstractVector{<:Real},
    points::AbstractVector{<:Real},
    partialelem::AbstractDict{<:AbstractVector{Int8},<:Real},
    ham::AbstractDict{<:AbstractVector{Int8},<:Real},
    generators::AbstractMatrix{Int8};
    atol::Real)::Dict{Vector{Int8},Float64}

    task = Threads.@spawn ham
    for i in eachindex(angles)
        partialelem = conjugate(partialelem, view(generators, :, i:i), -angles[i:i], atol=atol)
        partialham = fetch(task)
        _cost!(points, view(generators, :, i:i), partialelem, partialham, atol=atol)
        angles[i] = _minanglefind(points)
        task = Threads.@spawn conjugate(partialham, view(generators, :, i:i), -angles[i:i], atol=atol)
        # cosines[i] = cos(2 * angles[i])
        # sines[i] = -sin(2 * angles[i])
    end
    return fetch(task)
end

function errorfind!(ham::AbstractDict{<:AbstractVector{Int8},<:Real}, subalgebra::AbstractMatrix{Int8})::Float64
    """
    errorfind
    ----------
    Finds the error of the Hamiltonian after the transformation.
    """
    errornorm = 0.0
    fullnorm = 0.0
    for (key, value) in ham
        fullnorm += abs2(value)
        for string in eachcol(subalgebra)
            if !pauliprod(key, string)[3]
                errornorm += abs2(pop!(ham, key))
                break
            end
        end
    end
    return errornorm / fullnorm
end

function optimizer(ham::AbstractDict{<:AbstractVector{Int8},<:Real},
    subalgebra::AbstractMatrix{Int8},
    generators::AbstractMatrix{Int8},
    initangles::AbstractVector{<:Real}=pi * rand(size(generators, 2));
    method::AbstractString="roto",
    maxiter::Integer=0,
    convergence_tol::Real=1e-6,
    coeff_tol::Real=0,
    toltype::AbstractString="relerror",
    itertrack::Bool=false,
    timetrack::Bool=false)

    length(initangles) == size(generators, 2) || throw(ArgumentError("Incorrect number of initial angles. Expected $(size(generators, 2)), got $(length(initangles))."))

    irr = _mutirr(size(subalgebra, 2))
    subalgelem = Dict{Vector{Int8},Float64}((subalgebra[:, i], irr[i]) for i in axes(subalgebra, 2))
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
            println("Begin iteration $iter...")
            partialelem = conjugate(subalgelem, generators, angles, atol=coeff_tol)
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
                        angles = pi * rand(size(generators, 2))
                        t = time()
                        errorcache = 1.0
                        # angles +=  convergence_tol * randn(size(angles))

                    else
                        errorcache = relerror
                        println("Relative error after $iter iterations: $(sqrt(relerror))")
                        println()
                    end
                end
            end
        end
    end
end