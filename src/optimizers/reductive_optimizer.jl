function errorfind(ham::PauliSentence, subalgebra::PauliList)::Float64
    errornorm = 0.0
    fullnorm = 0.0
    for (key, value) in ham
        fullnorm += abs2(value)
        for h in subalgebra
            com(h, key, ham.qubits) == 0 || (errornorm += abs2(value); break)
        end
    end
    return errornorm / fullnorm
end

function _filter(ham::PauliSentence, subalgebra::PauliList)
    filteredham = copy(ham)
    for key in keys(ham)
        for h in subalgebra
            if com(h, key, ham.qubits) != 0
                delete!(filteredham, key)
                break
            end
        end
    end
    return filteredham
end

function _reductive_optimizer_step(
    ham::PauliSentence,
    abstrings::PauliList,
    generators::PauliList,
    angles::AbstractVector{<:Real},
    previous_ham::PauliSentence;
    method::Symbol,
    maxiter::Integer,
    convergence_tol::Real,
    coeff_tol::Real,
    toltype::Symbol,
    itertrack::Bool,
    timetrack::Bool,
)

    length(angles) == length(generators) || throw(
        ArgumentError(
            "Incorrect number of initial angles. Expected $(length(generators)),
            got $(length(angles)).",
        ),
    )

    subalgelem = PauliSentence(
        Dict{keytype(ham),Float64}(abstrings[end] => im^county(abstrings[end], ham.qubits)),
        ham.qubits,
    )
    if method == :roto
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
                    transformedham = ad(
                        previous_ham,
                        reverse(generators),
                        -reverse(angles),
                        atol=coeff_tol,
                    )
                    relerror = errorfind(transformedham, abstrings)
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
                                    :iterations => iter,
                                    :calls => 3 * iter * length(angles),
                                    :time => time() - t,
                                )
                            else
                                return Dict(
                                    :H => transformedham,
                                    :angles => angles,
                                    :iterations => iter,
                                    :calls => 3 * iter * length(angles),
                                )
                            end
                        else
                            if timetrack
                                return Dict(
                                    :H => transformedham,
                                    :angles => angles,
                                    :time => time() - t,
                                )
                            else
                                return Dict(:H => transformedham, :angles => angles)
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

function reductive_optimizer(
    ham::PauliSentence,
    abstrings::PauliList,
    symgenerators::AbstractVector{<:PauliList},
    initangles::AbstractVector{<:AbstractVector{<:Real}}=[
        pi * rand(length(symgen)) for symgen in symgenerators
    ];
    method::Symbol=:roto,
    maxiter::Integer=0,
    convergence_tol::Real=1e-6,
    coeff_tol::Real=0,
    toltype::Symbol=:relerror,
    itertrack::Bool=false,
    timetrack::Bool=false,
)

    angles = copy(initangles)
    relerror = 1.0
    iter = 0
    calls = 0
    stepham = ham
    previous_ham = ham
    t = 0.0
    for i in eachindex(abstrings)
        # ittol = max(1e-4, convergence_tol * 10.0^(1 - i))
        println("Begin optimization for abelian element $i")
        opt = _reductive_optimizer_step(
            stepham,
            abstrings[1:i],
            symgenerators[i],
            angles[i],
            previous_ham,
            method=method,
            maxiter=maxiter,
            convergence_tol=convergence_tol,
            coeff_tol=coeff_tol,
            toltype=toltype,
            itertrack=itertrack,
            timetrack=timetrack,
        )
        println()
        previous_ham = opt[:H]
        stepham = _filter(opt[:H], abstrings[1:i])
        angles[i] = opt[:angles]
        if itertrack
            iter += opt[:iterations]
            calls += opt[:calls]
        end
        timetrack && (t += opt[:time])
    end
    finalham = ad(
        ham,
        reverse!(vcat(symgenerators...)),
        -reverse!(vcat(angles...)),
        atol=coeff_tol,
    )
    relerror = errorfind!(finalham, abstrings)
    println("Combined relative error: $(sqrt(relerror))")

    println("Optimization complete. Combined relative error: $(sqrt(relerror))")
    if itertrack
        if timetrack
            return Dict(
                :H => finalham,
                :angles => angles,
                :error => sqrt(relerror),
                :iterations => iter,
                :calls => calls,
                :time => t,
            )
        else
            return Dict(
                :H => finalham,
                :angles => angles,
                :error => sqrt(relerror),
                :iterations => iter,
                :calls => calls,
            )
        end
    else
        if timetrack
            return Dict(
                :H => finalham,
                :angles => angles,
                :error => sqrt(relerror),
                :time => t,
            )
        else
            return Dict(:H => finalham, :angles => angles, :error => sqrt(relerror))
        end
    end
end