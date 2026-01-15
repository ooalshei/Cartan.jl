function errorfind(ham::PauliSentence, subalgebra::PauliList)::Float64
    errornorm = 0.0
    fullnorm = 0.0
    for (key, value) in ham
        fullnorm += abs2(value)
        for h in subalgebra
            com(h, key, ham.qubits).second || (errornorm += abs2(value); break)
        end
    end
    return errornorm / fullnorm
end

function _filter(ham::PauliSentence, subalgebra::PauliList)
    filteredham = copy(ham)
    for key in keys(ham)
        for h in subalgebra
            if !com(h, key, ham.qubits).second
                delete!(filteredham, key)
                break
            end
        end
    end
    return filteredham
end

function _reductive_optimizer_step(
    ham::PauliSentence,
    bstrings::PauliList,
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
        Dict{keytype(ham),Float64}(bstrings[end] => im^county(bstrings[end], ham.qubits)),
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
                    relerror = errorfind(transformedham, bstrings)
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

@doc raw"""
    reductive_optimizer

Perform a reductive KHK optimization.

# Arguments
- `ham::PauliSentence`: the Hamiltonian to be transformed.
- `bstrings::PauliList`: the strings that generate (under multiplication) the target Abelian
subalgebra to rotate into.
- `symgenerators::AbstractVector{<:PauliList}`: the elements of the ``\mathfrak{k}``,
fragmented with respect to `bstrings`. See [`fragmentedsubspaces`](@ref).
- `initangles::AbstractVector{<:AbstractVector{<:Real}}=[
    pi * rand(length(symgen)) for symgen in symgenerators
]`: initial angles for the generators.
- `method::Symbol=:roto`: optimization method to use (currently only `:roto` is supported).
- `maxiter::Integer=0`: maximum number of iterations (0 for unlimited).
- `convergence_tol::Real=1e-6`: tolerance for convergence.
- `coeff_tol::Real=0`: tolerance for ignoring small coefficients in the Pauli sentences.
- `toltype::Symbol=:relerror`: type of tolerance to use (currently only `:relerror` is
supported).
- `itertrack::Bool=false`: whether to track the number of iterations and calls.
- `timetrack::Bool=false`: whether to track the total time taken.

Rotate `ham` reductively into the Cartan subalgebra using the fragmented generators in
`symgenerators`. `ham` is first transformed so that it commutes with the first element in
`bstrings`. Then it is subsequently transformed until it commutes with all elements in
`bstrings`. The end result lies in the Cartan subalgebra. The optimization starts from
`initangles` and proceeds until the error is below `convergence_tol` or `maxiter` is
reached. Small coefficients below `coeff_tol` are ignored during the optimization. If
`itertrack` is true, the number of iterations and function calls are tracked. If `timetrack`
is true, the total time taken is tracked.

# Notes
- Currently, only the [rotosolve](https://quantum-journal.org/papers/q-2021-01-28-391/)
optimizer is supported. This optimizer iteratively updates the angles for each generator
individually to minimize the cost function. This method is crucial for the quantum-assisted
optimization algorithm described in the [paper](https://arxiv.org/abs/2512.06070).
- Currently, the only available option for `toltype` is `:relerror`, which uses the relative
error between `ham` and its projection onto the Cartan `subalgebra` as the convergence
metric. This is the ratio of the Hilbert-Schmidt norm of the part of the transformed `ham`
that lies outside of the `subalgebra` to the Hilbert-Schmidt norm of `ham`.

See also [`optimizer`](@ref).
"""
function reductive_optimizer(
    ham::PauliSentence,
    bstrings::PauliList,
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
    for i in eachindex(bstrings)
        # ittol = max(1e-4, convergence_tol * 10.0^(1 - i))
        println("Begin optimization for abelian element $i")
        opt = _reductive_optimizer_step(
            stepham,
            bstrings[1:i],
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
        stepham = _filter(opt[:H], bstrings[1:i])
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
    relerror = errorfind!(finalham, bstrings)
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