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
function _mutirr(len::Int; seed::Irrational=pi)::Vector{Float64}
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

function _minanglefind(points::Vector{Float64})::Float64
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
    return 0.5 * atan(points[2] - (points[1] + points[3]) / 2, (points[1] - points[3]) / 2) + pi / 2

end

function cost(generators::Matrix{Int8},
    angles::Vector{Float64},
    subalgelem::Dict{Vector{Int8},Float64},
    ham::Dict{Vector{Int8},Float64};
    tol::Float64=0.0)::Float64
    raw"""
    Calculates the cost function :math:`\mathrm{Tr}(KvK^\dag, \mathcal{H})`.
    """
    sentence = conjugate(subalgelem, generators, angles, tol=tol)
    costval = 0.0
    for (key, value) in ham
        key in keys(sentence) && (costval += value * sentence[key])
    end
    return costval
end

function _cost!(ret::Vector{Float64},
    generators::Matrix{Int8},
    angles::Vector{Float64},
    partialelem::Dict{Vector{Int8},Float64},
    ham::Dict{Vector{Int8},Float64};
    tol::Float64=0.0)::Nothing
    raw"""
    Calculates the cost function :math:`\mathrm{Tr}(KvK^\dag, \mathcal{H})` at 0, pi/4, pi/2.
    """

    angles[end] = 0
    ret[1] = cost(generators, angles, partialelem, ham, tol=tol)
    angles[end] = pi / 4
    ret[2] = cost(generators, angles, partialelem, ham, tol=tol)
    angles[end] = pi / 2
    ret[3] = cost(generators, angles, partialelem, ham, tol=tol)
    nothing
end

function _rotostep!(partialelem::Dict{Vector{Int8},Float64},
    points::Vector{Float64},
    angles::Vector{Float64},
    generators::Matrix{Int8},
    ham::Dict{Vector{Int8},Float64};
    tol::Float64)::Nothing

    for i in eachindex(angles)
        partialelem = _conjugate(partialelem, generators[:, i], -angles[i], tol=tol)
        _cost!(points, generators[:, 1:i], angles[1:i], partialelem, ham, tol=tol)
        angles[i] = _minanglefind(points)
        # cosines[i] = cos(2 * angles[i])
        # sines[i] = -sin(2 * angles[i])
    end
    nothing
end

function errorfind!(ham::Dict{Vector{Int8},Float64}, subalgebra::Matrix{Int8})::Float64
    """
    errorfind
    ----------
    Finds the error of the Hamiltonian after the transformation.
    """
    errornorm = 0.0
    fullnorm = 0.0
    for key in keys(ham)
        fullnorm += abs2(ham[key])
        for i in axes(subalgebra, 2)
            if !pauliprod(key, subalgebra[:, i])[3]
                errornorm += abs2(pop!(ham, key))
                break
            end
        end
    end
    return errornorm / fullnorm
end

function optimizer(ham::Dict{Vector{Int8},Float64},
    subalgebra::Matrix{Int8},
    generators::Matrix{Int8},
    initangles::Vector{Float64}=pi * rand(size(generators, 2));
    method::String="roto",
    maxiter::Int=0,
    mintol::Float64=1e-6,
    tol::Float64=0.0,
    toltype::String="relerror")::Dict{String,Union{Dict{Vector{Int8},Float64},Vector{Float64}}}

    length(initangles) == size(generators, 2) || ArgumentError("Incorrect number of initial angles. Expected $(size(generators, 2)), got $(length(initangles)).")

    irr = _mutirr(size(subalgebra, 2))
    subalgelem = Dict{Vector{Int8},Float64}((subalgebra[:, i], irr[i]) for i in axes(subalgebra, 2))
    if method == "roto"
        angles = copy(initangles)
        # cosines = cos.(2 .* angles)
        # sines = .-sin.(2 .* angles)
        points = zeros(Float64, 3)
        errorcache = 1.0

        iter = 0
        while true
            iter += 1
            println("Begin iteration $iter...")
            partialelem = conjugate(subalgelem, generators, angles, tol=tol)
            _rotostep!(partialelem, points, angles, generators, ham, tol=tol)
            if (iter % 10 == 0) | (iter == maxiter)
                # println(cost(generators, angles, subalgelem, ham, tol=tol))
                if toltype == "relerror"
                    transformedham = conjugate(ham, reverse(generators, dims=2), -reverse(angles), tol=tol)
                    relerror = errorfind!(transformedham, subalgebra)
                    if (relerror <= mintol^2) | (iter == maxiter)
                        iter == maxiter ? println("Max iterations reached.") : println("Converged in $iter iterations.")
                        println("Final relative error: $(sqrt(relerror))")
                        return Dict("H" => transformedham, "angles" => angles)
                    elseif abs(sqrt(errorcache) - sqrt(relerror)) < mintol / 10
                        println("Relative error after $iter iterations: $(sqrt(relerror))")
                        println("Convergence too slow. Starting over with another initial guess.")
                        println()
                        iter = 0
                        errorcache = 1.0
                        angles = pi * rand(size(generators, 2))

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
