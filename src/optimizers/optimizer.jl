function _mutirr(len::Integer; seed::Irrational=pi)::Vector{Float64}
    result = Vector{Float64}(undef, len)
    result[1] = seed % 1
    for i in 2:len
        result[i] = (seed * result[i-1]) % 1
    end
    return result
end

_minanglefind(points::AbstractVector{<:Real})::Float64 =
    0.5 * atan(2 * points[2] - (points[1] + points[3]), (points[1] - points[3])) + pi / 2

function cost(
    subalgelem::PauliSentence,
    generators::Union{UPauli,PauliList},
    angles::Union{Real,AbstractVector{<:Real}},
    ham::PauliSentence;
    atol::Real=0.0,
)::Float64

    sentence = ad(subalgelem, generators, angles, atol=atol)
    result::Float64 = 0.0
    for (key, value) in ham
        haskey(sentence, key) &&
            (result += (-1)^county(key, ham.qubits) * value * sentence[key])
    end
    return result
end

@inline function _cost!(
    ret::AbstractVector{<:Real},
    generator::UPauli,
    partialelem::PauliSentence,
    partialham::PauliSentence;
    atol::Real,
)

    @inbounds Threads.@threads for i in 0:2
        ret[i+1] = cost(partialelem, generator, i * pi / 4, partialham, atol=atol)
    end
    nothing
end

function _rotostep!(
    angles::AbstractVector{<:Real},
    points::AbstractVector{<:Real},
    partialelem::PauliSentence,
    ham::PauliSentence,
    generators::PauliList;
    atol::Real,
)

    task = ham
    @inbounds for i in eachindex(angles)
        pauligen = UPauli(generators[i], ham.qubits)
        ad!(partialelem, pauligen, -angles[i], atol=atol)
        partialham = fetch(task)
        _cost!(points, pauligen, partialelem, partialham, atol=atol)
        angles[i] = _minanglefind(points)
        task = Threads.@spawn ad!(partialham, pauligen, -angles[i], atol=atol)
    end
    return fetch(task)
end

function errorfind!(ham::PauliSentence, subalgebra::PauliList)::Float64
    errornorm = 0.0
    fullnorm = 0.0
    for (key, value) in ham
        fullnorm += abs2(value)
        # unique!(com.(subalgebra, key, ham.qubits)) == [0] ||
        #     (errornorm += abs2(value); delete!(ham, key))
        for h in subalgebra
            if com(h, key, ham.qubits) != 0
                errornorm += abs2(value)
                delete!(ham, key)
                break
            end
        end
    end
    return errornorm / fullnorm
end

function optimizer(
    ham::PauliSentence,
    subalgebra::PauliList,
    generators::PauliList,
    initangles::AbstractVector{<:Real}=pi * rand(length(generators));
    method::Symbol=:roto,
    maxiter::Integer=0,
    convergence_tol::Real=1e-6,
    coeff_tol::Real=0,
    toltype::Symbol=:relerror,
    itertrack::Bool=false,
    timetrack::Bool=false,
)

    length(initangles) == length(generators) || throw(
        ArgumentError(
            "Incorrect number of initial angles. Expected $(length(generators)),
            got $(length(initangles)).",
        ),
    )

    irr = _mutirr(length(subalgebra))
    subalgelem = PauliSentence(subalgebra, (im) .^ county.(subalgebra, ham.qubits) .* irr)
    if method == :roto
        angles = copy(initangles)
        points = Vector{Float64}(undef, 3)
        # errorcache = 1.0

        iter = 0
        t = time()
        while true
            iter += 1
            # println("Begin iteration $iter...")
            partialelem = ad(subalgelem, generators, angles, atol=coeff_tol)
            transformedham = PauliSentence{keytype(ham),ComplexF64}(ham)
            _rotostep!(
                angles,
                points,
                partialelem,
                transformedham,
                generators,
                atol=coeff_tol,
            )
            if (iter % 10 == 0) | (iter == maxiter)
                if toltype == :relerror
                    relerror = errorfind!(transformedham, subalgebra)
                    if (relerror <= convergence_tol^2) | (iter == maxiter)
                        if iter == maxiter
                            println("Max iterations reached.")
                        else
                            println("Converged in $iter iterations.")
                        end
                        println("Final relative error: $(sqrt(relerror))")

                        if itertrack
                            if timetrack
                                return Dict(
                                    :H => transformedham,
                                    :angles => angles,
                                    :error => sqrt(relerror),
                                    :iterations => iter,
                                    :calls => 3 * iter * length(angles),
                                    :time => time() - t,
                                )
                            else
                                return Dict(
                                    :H => transformedham,
                                    :angles => angles,
                                    :error => sqrt(relerror),
                                    :iterations => iter,
                                    :calls => 3 * iter * length(angles),
                                )
                            end

                        else
                            if timetrack
                                return Dict(
                                    :H => transformedham,
                                    :angles => angles,
                                    :error => sqrt(relerror),
                                    :time => time() - t,
                                )
                            else
                                return Dict(
                                    :H => transformedham,
                                    :angles => angles,
                                    :error => sqrt(relerror),
                                )
                            end
                        end

                    # elseif (sqrt(errorcache) - sqrt(relerror)) <= 0
                    #     println("Relative error after $iter iterations: $(sqrt(relerror))")
                    #     println("Convergence too slow. Starting over.")
                    #     println()
                    #     iter = 0
                    #     angles = pi * rand(length(generators))
                    #     t = time()
                    #     errorcache = 1.0
                    # angles +=  convergence_tol * randn(size(angles))

                    else
                        # errorcache = relerror
                        iter % 50 == 0 && println(
                            "Relative error after $iter iterations: $(sqrt(relerror))",
                        )
                        println()
                    end
                end
            end
        end
    end
end