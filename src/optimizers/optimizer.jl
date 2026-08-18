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

@doc raw"""
    _costcoeffs(subalgelem, generator, ham)

Coefficients ``(c_1, c_2)`` of the single-generator cost
``\langle H, e^{-i\theta g} A e^{i\theta g}\rangle = c_0 + c_1\cos 2\theta + c_2\sin 2\theta``
for the Cartan element `A = subalgelem` and the generator `g = generator`.

Rotating about one generator only mixes each Pauli `P` with its partner `g ⊻ P`, so both
coefficients follow from a single sweep over `ham`; the constant `c_0` never enters
[`_minanglefind`](@ref) and is not accumulated.
"""
function _costcoeffs(
    subalgelem::PauliSentence,
    generator::UPauli{<:Unsigned,Q},
    ham::PauliSentence{<:Unsigned,<:Number,Q},
) where {Q}

    g = generator.string
    modsine = _ipow(county(g, Q) + 1)
    c1 = zero(ComplexF64)
    c2 = zero(ComplexF64)
    for (key, value) in ham
        product = com(g, key, Q)
        product.second && continue  # untouched by the rotation: contributes to c₀ only
        weight = isodd(county(key, Q)) ? -value : value
        c1 += weight * get(subalgelem, key, zero(ComplexF64))
        c2 +=
            weight *
            modsine *
            _prodsign(g, product.first, Q) *
            get(subalgelem, product.first, zero(ComplexF64))
    end
    return real(c1), real(c2)
end

# Minimum of c₀ + c₁cos(2θ) + c₂sin(2θ). Equivalent to `_minanglefind` on the three samples
# θ = 0, π/4, π/2, which are exactly (c₀+c₁, c₀+c₂, c₀-c₁).
_minanglefind(c1::Real, c2::Real)::Float64 = 0.5 * atan(c2, c1) + pi / 2

function _serial_rotostep!(
    angles::AbstractVector{<:Real},
    partialelem::PauliSentence,
    ham::PauliSentence{<:Unsigned,<:Number,Q},
    generators::PauliList{T,Q};
    atol::Real,
) where {T,Q}

    partialham = ham
    @inbounds for i in eachindex(angles)
        pauligen = UPauli{T,Q}(generators[i])
        ad!(partialelem, pauligen, -angles[i], atol=atol)
        angles[i] = _minanglefind(_costcoeffs(partialelem, pauligen, partialham)...)
        ad!(partialham, pauligen, -angles[i], atol=atol)
    end
    return partialham
end

function _parallel_rotostep!(
    angles::AbstractVector{<:Real},
    partialelem::PauliSentence,
    ham::PauliSentence{<:Unsigned,<:Number,Q},
    generators::PauliList{T,Q};
    atol::Real,
) where {T,Q}

    task = ham
    @inbounds for i in eachindex(angles)
        pauligen = UPauli{T,Q}(generators[i])
        ad!(partialelem, pauligen, -angles[i], atol=atol)
        partialham = fetch(task)
        angles[i] = _minanglefind(_costcoeffs(partialelem, pauligen, partialham)...)
        task = Threads.@spawn ad!(partialham, pauligen, -angles[i], atol=atol)
    end
    return fetch(task)
end

# Resolved on every call: the number of threads is a property of the running session, not of
# the session that precompiled the package.
_rotostep!(args...; kwargs...) = if isone(Threads.nthreads())
    _serial_rotostep!(args...; kwargs...)
else
    _parallel_rotostep!(args...; kwargs...)
end

"""
    _commutes(key, subalgebra, Q)

Whether the Pauli `key` commutes with every element of `subalgebra`, i.e. whether it lies in
the centralizer the optimization rotates the Hamiltonian into.
"""
_commutes(key::Unsigned, subalgebra::PauliList, Q::Integer) =
    all(h -> com(h, key, Q).second, subalgebra)

"""
    errorfind!(ham::PauliSentence, subalgebra::PauliList)

Delete from `ham` every term that does not commute with all of `subalgebra`, and return the
squared relative Hilbert-Schmidt weight of what was deleted.
"""
function errorfind!(
    ham::PauliSentence{<:Unsigned,<:Number,Q},
    subalgebra::PauliList,
)::Float64 where {Q}
    errornorm = 0.0
    fullnorm = 0.0
    # `filter!` deletes in one pass and, unlike deleting while iterating, is defined
    # behaviour for the underlying dictionary.
    filter!(ham) do (key, value)
        weight = abs2(value)
        fullnorm += weight
        _commutes(key, subalgebra, Q) || (errornorm += weight; return false)
        return true
    end
    return errornorm / fullnorm
end

@doc raw"""
    optimizer

Perform the KHK optimization.

# Arguments
- `ham::PauliSentence`: the Hamiltonian to be transformed.
- `subalgebra::PauliList`: the target Abelian subalgebra to rotate into.
- `generators::PauliList`: the elements of the ``\mathfrak{k}`` subalgebra.
- `initangles::AbstractVector{<:Real}=pi * rand(length(generators))`: initial angles for the
generators.
- `method::Symbol=:roto`: optimization method to use (currently only `:roto` is supported).
- `maxiter::Integer=0`: maximum number of iterations (0 for unlimited).
- `convergence_tol::Real=1e-6`: tolerance for convergence.
- `coeff_tol::Real=0`: tolerance for ignoring small coefficients in the Pauli sentences.
- `toltype::Symbol=:relerror`: type of tolerance to use (currently only `:relerror` is
supported).
- `itertrack::Bool=false`: whether to track the number of iterations and calls.
- `timetrack::Bool=false`: whether to track the total time taken.

Rotate `ham` into the specificed Cartan `subalgebra` using `generators`. The optimization
starts from `initangles` and proceeds until the error is below `convergence_tol` or
`maxiter` is reached. Small coefficients below `coeff_tol` are ignored during the
optimization. If `itertrack` is true, the number of iterations and function calls are
tracked. If `timetrack` is true, the total time taken is tracked.

# Notes
- Currently, only the [rotosolve](https://quantum-journal.org/papers/q-2021-01-28-391/)
optimizer is supported. This optimizer iteratively updates the angles for each generator
individually to minimize the cost function.
- Currently, the only available option for `toltype` is `:relerror`, which uses the relative
error between `ham` and its projection onto the Cartan `subalgebra` as the convergence
metric. This is the ratio of the Hilbert-Schmidt norm of the part of the transformed `ham`
that lies outside of the `subalgebra` to the Hilbert-Schmidt norm of `ham`.

See also [`reductive_optimizer`](@ref).
"""
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
        ArgumentError("Incorrect number of initial angles. Expected $(length(generators)),
                      got $(length(initangles))."),
    )

    irr = _mutirr(length(subalgebra))
    subalgelem = PauliSentence(subalgebra, (im) .^ county.(subalgebra, ham.qubits) .* irr)
    # Converted once: every sweep starts from a fresh copy of the Hamiltonian, and copying a
    # dictionary that already holds the right element type is a bulk copy rather than a
    # rehash of every term.
    complexham = PauliSentence{keytype(ham),ComplexF64}(ham)
    if method == :roto
        angles = copy(initangles)
        # errorcache = 1.0

        iter = 0
        t = time()
        while true
            iter += 1
            # println("Begin iteration $iter...")
            partialelem = ad(subalgelem, generators, angles, atol=coeff_tol)
            transformedham = copy(complexham)
            _rotostep!(angles, partialelem, transformedham, generators, atol=coeff_tol)
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
                                return (
                                    H=transformedham,
                                    generators=generators,
                                    angles=angles,
                                    error=sqrt(relerror),
                                    iterations=iter,
                                    calls=3 * iter * length(angles),
                                    time=time() - t,
                                )
                            else
                                return (
                                    H=transformedham,
                                    generators=generators,
                                    angles=angles,
                                    error=sqrt(relerror),
                                    iterations=iter,
                                    calls=3 * iter * length(angles),
                                )
                            end

                        else
                            if timetrack
                                return (
                                    H=transformedham,
                                    generators=generators,
                                    angles=angles,
                                    error=sqrt(relerror),
                                    time=time() - t,
                                )
                            else
                                return (
                                    H=transformedham,
                                    generators=generators,
                                    angles=angles,
                                    error=sqrt(relerror),
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
