function generatex(n::Integer, ::Type{T}=UInt) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n))
    number = T(1)
    for i in 1:n
        result[i] = number
        number <<= 1
    end
    return result
end

function generatez(n::Integer, ::Type{T}=UInt) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n))
    number = T(1) << n
    for i in 1:n
        result[i] = number
        number <<= 1
    end
    return result
end

function generatey(n::Integer, ::Type{T}=UInt) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n))
    number = (T(1) << n) | T(1)
    for i in 1:n
        result[i] = number
        number <<= 1
    end
    return result
end

function generatexx(n::Integer, ::Type{T}=UInt; pbc::Bool=false) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n - 1))
    number = T(3)
    for i in 1:(n-1)
        result[i] = number
        number <<= 1
    end
    pbc && return push!(result, (T(1) << (n - 1)) | T(1))
    return result
end

function generatezz(n::Integer, ::Type{T}=UInt; pbc::Bool=false) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n - 1))
    number = (T(1) << (n + 1)) | (T(1) << n)
    for i in 1:(n-1)
        result[i] = number
        number <<= 1
    end
    pbc && return push!(result, (T(1) << (2 * n - 1)) | (T(1) << n))
    return result
end

function generateyy(n::Integer, ::Type{T}=UInt; pbc::Bool=false) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n - 1))
    number = (T(1) << (n + 1)) | (T(1) << n) | T(3)
    for i in 1:(n-1)
        result[i] = number
        number <<= 1
    end
    pbc &&
        return push!(result, (T(1) << (2 * n - 1)) | (T(1) << n) | (T(1) << (n - 1)) | T(1))
    return result
end

function generatexy(n::Integer, ::Type{T}=UInt; pbc::Bool=false) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n - 1))
    number = (T(1) << n) | T(3)
    for i in 1:(n-1)
        result[i] = number
        number <<= 1
    end
    pbc && return push!(result, (T(1) << (2 * n - 1)) | (T(1) << (n - 1)) | T(1))
    return result
end

function generateyx(n::Integer, ::Type{T}=UInt; pbc::Bool=false) where {T<:Unsigned}
    result = PauliList{T,n}(zeros(T, n - 1))
    number = (T(1) << (n + 1)) | T(3)
    for i in 1:(n-1)
        result[i] = number
        number <<= 1
    end
    pbc && return push!(result, (T(1) << n) | (T(1) << (n - 1)) | T(1))
    return result
end

"""
    hamiltonian(model, n, couplings, T=UInt; nf=0, pbc=false)

Construct common lattice Hamiltonians.

Output a `PauliSentence` representing the Hamiltonian specified by `model` on `n` qubits
with given `couplings`. Periodic boundary conditions can be included with `pbc=true`.
The `"GN"` and `"SYK"` models support an additional `nf` argument specifying the number of
flavors.
"""
function hamiltonian(
    model::AbstractString,
    n::Integer,
    couplings::AbstractVector{<:Real}=[1.0],
    ::Type{T}=UInt;
    nf::Integer=0,
    pbc::Bool=false,
) where {T<:Unsigned}
    if uppercase(model) == "ISING"
        length(couplings) == 1 || throw(
            ArgumentError(
                "Incorrect number of couplings. Expected 1 (J), got $(length(couplings)).
                H = -JXX",
            ),
        )
        strings = generatexx(n, T, pbc=pbc)
        coefficients = fill(-couplings[1], length(strings))

    elseif uppercase(model) == "XY"
        length(couplings) == 1 || throw(
            ArgumentError(
                "Incorrect number of couplings. Expected 1 (J), got $(length(couplings)).
                H = -J(XX + YY)",
            ),
        )
        strings = append!(generatexx(n, T, pbc=pbc), generateyy(n, T, pbc=pbc))
        coefficients = fill(-couplings[1], length(strings))

    elseif uppercase(model) == "TFIM"
        if length(couplings) != 2
            throw(ArgumentError("Incorrect number of couplings. Expected 2 (J,g),
                                got $(length(couplings)).H = -J(XX + gZ)"))
        else
            strings = append!(generatexx(n, T, pbc=pbc), generatez(n, T))
            coefficients = fill(-couplings[1], length(strings))
            coefficients[(end-n+1):end] *= couplings[2]
        end

    elseif uppercase(model) == "TFXY"
        if length(couplings) != 2
            throw(ArgumentError("Incorrect number of couplings. Expected 2 (J,g),
                                got $(length(couplings)). H = -J(XX + YY + gZ)"))
        else
            strings = append!(
                generatexx(n, T, pbc=pbc),
                generateyy(n, T, pbc=pbc),
                generatez(n, T),
            )
            coefficients = fill(-couplings[1], length(strings))
            coefficients[(end-n+1):end] *= couplings[2]
        end

    elseif uppercase(model) == "MFIM"
        if length(couplings) != 3
            throw(ArgumentError("Incorrect number of couplings. Expected 3 (J,gh,gl),
                                got $(length(couplings)). H = -J(XX + gh Z + gl X)"))
            # elseif couplings[2] == 0
            #     strings, coefficients = hamiltonian("TFIM", n, couplings[1:2], pbc=pbc)
        else
            strings = append!(generatexx(n, T, pbc=pbc), generatez(n, T), generatex(n, T))
            coefficients = fill(-couplings[1], length(strings))
            coefficients[(end-2*n+1):(end-n)] *= couplings[2]
            coefficients[(end-n+1):end] *= couplings[3]
        end

    elseif uppercase(model) == "HIM"
        if length(couplings) != 4
            throw(ArgumentError("Incorrect number of couplings. Expected 4 (J,gh,gl,ℓ),
                                got $(length(couplings)). H = -J(ℓ)(XX + gh Z + gl X)"))
        else
            strings = append!(generatexx(n, T, pbc=pbc), generatez(n, T), generatex(n, T))
            coefficients = fill(-couplings[1], length(strings))
            η = [cosh(-couplings[4] + i * 2 * couplings[4] / (n - 1)) for i in 1:n]  # η_i
            x = pbc ? n : n - 1
            for i in 1:(n-1)
                coefficients[i] *= (η[i] + η[i+1] / 2)  # XX couplings
                coefficients[i+x] *= couplings[2] * η[i]  # Z couplings
                coefficients[i+x+n] *= couplings[3] * η[i]  # X couplings
            end
            pbc && (coefficients[n] *= (η[1] + η[end]) / 2)  # XX coupling
            coefficients[end-n] *= couplings[2] * η[end]  # Z coupling
            coefficients[end] *= couplings[3] * η[end]  # X coupling
        end

    elseif uppercase(model) == "HEISENBERG"
        length(couplings) == 1 ||
            throw(ArgumentError("Incorrect number of couplings. Expected 1 (J),
                                got $(length(couplings)). H = -J(XX + YY + ZZ)"))
        strings = append!(
            generatexx(n, T, pbc=pbc),
            generateyy(n, T, pbc=pbc),
            generatezz(n, T, pbc=pbc),
        )
        coefficients = fill(-couplings[1], length(strings))

    elseif uppercase(model) == "XXZ"
        length(couplings) == 2 ||
            throw(ArgumentError("Incorrect number of couplings. Expected 2 (J,Δ),
                                got $(length(couplings)). H = -J(XX + YY + Δ ZZ)"))
        strings = append!(
            generatexx(n, T, pbc=pbc),
            generateyy(n, T, pbc=pbc),
            generatezz(n, T, pbc=pbc),
        )
        coefficients = fill(-couplings[1], length(strings))
        coefficients[(2*length(strings)÷3+1):end] *= couplings[2]

    elseif uppercase(model) == "XXZ_EXT"
        length(couplings) == 3 ||
            throw(ArgumentError("Incorrect number of couplings. Expected 3 (J,Δ,g),
                                got $(length(couplings)). H = -J(XX + YY + Δ ZZ) + gZ"))
        strings = append!(
            generatexx(n, T, pbc=pbc),
            generateyy(n, T, pbc=pbc),
            generatezz(n, T, pbc=pbc),
            generatez(n, T),
        )
        coefficients = fill(-couplings[1], length(strings))
        coefficients[(2*(length(strings)-n)÷3+1):(end-n)] *= couplings[2]
        coefficients[(end-n+1):end] .= couplings[3]

    elseif uppercase(model) == "GN"
        length(couplings) == 3 || throw(
            ArgumentError(
                "Incorrect number of couplings. Expected 3 (μ,G,m),
                got $(length(couplings)). H = (1-2μ)(YX - XY) + (2m-4nf nG)Z - 4G ZZ ",
            ),
        )
        nf > 0 || throw(ArgumentError("Number of flavors must be greater than 0."))
        strings = generateyx(nf * n, T)
        setdiff!(
            strings,
            T[
                (T(1) << (i * n - 1)) | (T(1) << (i * n)) | (T(1) << (i * n + n * nf)) for
                i in 1:(nf-1)
            ],
        )
        temp = generatexy(nf * n, T)
        setdiff!(
            temp,
            T[
                (T(1) << (i * n - 1)) | (T(1) << (i * n)) | (T(1) << (i * n + n * nf - 1))
                for i in 1:(nf-1)
            ],
        )
        append!(strings, temp)
        temp = generatez(nf * n, T)
        append!(strings, temp)
        temp .&= (T(1) << (nf * n + n)) - T(1)
        setdiff!(temp, [0])
        for i in 0:(nf-1)
            temp2 = temp .<< (i * n)
            for j in 1:(nf-i-1)
                append!(strings, temp2 .| (temp2 .<< (j * n)))
            end
        end
        coefficients = [
            fill(1 - 2 * couplings[1], nf * (n - 1))
            fill(-(1 - 2 * couplings[1]), nf * (n - 1))
            fill(2 * couplings[3] - 4 * couplings[2] * (nf - 1), n * nf)
            fill(-4 * couplings[2], n * nf * (nf - 1) ÷ 2)
        ]

    elseif uppercase(model) == "SYK"

    else
        throw(ArgumentError("Model not recognized."))
    end
    return PauliSentence(strings, coefficients)
end
