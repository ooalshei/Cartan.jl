"""
cartan
------
This module contains methods that build Pauli strings given a Hamiltonian, generate the dynamical Lie algebra, and find
Cartan decompositions. The identity along with the three Pauli matrices are referred to as (1, 2, 3, 4).
"""
function generatexx(n::Integer; pbc::Bool=false)::Matrix{Int8}
    """
    generatexx
    ----------
    This function generates the XX strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] .= Int8(2)
    end
    if pbc
        strings = [strings ones(Int8, n)]
        strings[[n, 1], n] .= Int8(2)
    end
    return strings
end

function generateyy(n::Integer; pbc::Bool=false)::Matrix{Int8}
    """
    generateyy
    ----------
    This function generates the YY strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] .= Int8(3)
    end
    if pbc
        strings = [strings ones(Int8, n)]
        strings[[n, 1], n] .= Int8(3)
    end
    return strings
end

function generatezz(n::Integer; pbc::Bool=false)::Matrix{Int8}
    """
    generatezz
    ----------
    This function generates the ZZ strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] .= Int8(4)
    end
    if pbc
        strings = [strings ones(Int8, n)]
        strings[[n, 1], n] .= Int8(4)
    end
    return strings
end

function generatexy(n::Integer)::Matrix{Int8}
    """
    generatexy
    ----------
    This function generates the XY strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] = Int8[2, 3]
    end
    return strings
end

function generateyx(n::Integer)::Matrix{Int8}
    """
    generateyx
    ----------
    This function generates the YX strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] = Int8[3, 2]
    end
    return strings
end

function generatex(n::Integer)::Matrix{Int8}
    """
    generatex
    ---------
    This function generates the X strings for a chain of length n.
    """
    strings = ones(Int8, n, n)
    strings[CartesianIndex.(1:n, 1:n)] .= Int8(2)
    return strings
end

function generatey(n::Integer)::Matrix{Int8}
    """
    generatey
    ---------
    This function generates the Y strings for a chain of length n.
    """
    strings = ones(Int8, n, n)
    strings[CartesianIndex.(1:n, 1:n)] .= Int8(3)
    return strings
end

function generatez(n::Integer)::Matrix{Int8}
    """
    generatez
    ---------
    This function generates the Z strings for a chain of length n.
    """
    strings = ones(Int8, n, n)
    strings[CartesianIndex.(1:n, 1:n)] .= Int8(4)
    return strings
end

function hamiltonian(model::AbstractString, n::Integer, couplings::AbstractVector{<:Real}=[1.0]; nf::Integer=0, pbc::Bool=false)::Tuple{Matrix{Int8},Vector{Float64}}
    """
    hamiltonian
    -----------
    This function generates the Hamiltonian for a given model.
    """
    if uppercase(model) == "ISING"
        length(couplings) == 1 || throw(ArgumentError("Incorrect number of couplings. Expected 1 (J), got $(length(couplings)). H = -JXX"))
        strings = generatexx(n, pbc=pbc)
        coefficients = fill(-couplings[1], size(strings, 2))

    elseif uppercase(model) == "XY"
        length(couplings) == 1 || throw(ArgumentError("Incorrect number of couplings. Expected 1 (J), got $(length(couplings)). H = -J(XX + YY)"))
        strings = [generatexx(n, pbc=pbc) generateyy(n, pbc=pbc)]
        coefficients = fill(-couplings[1], size(strings, 2))

    elseif uppercase(model) == "TFIM"
        if length(couplings) != 2
            throw(ArgumentError("Incorrect number of couplings. Expected 2 (J,g), got $(length(couplings)). H = -J(XX + gZ)"))
        elseif couplings[2] == 0
            strings, coefficients = hamiltonian("Ising", n, [couplings[1]], pbc=pbc)
        else
            strings = [generatexx(n, pbc=pbc) generatez(n)]
            coefficients = fill(-couplings[1], size(strings, 2))
            coefficients[end-n+1:end] *= couplings[2]
        end

    elseif uppercase(model) == "TFXY"
        if length(couplings) != 2
            throw(ArgumentError("Incorrect number of couplings. Expected 2 (J,g), got $(length(couplings)). H = -J(XX + YY + gZ)"))
        elseif couplings[2] == 0
            strings, coefficients = hamiltonian("XY", n, [couplings[1]], pbc=pbc)
        else
            strings = [generatexx(n, pbc=pbc) generateyy(n, pbc=pbc) generatez(n)]
            coefficients = fill(-couplings[1], size(strings, 2))
            coefficients[end-n+1:end] *= couplings[2]
        end

    elseif uppercase(model) == "CFIM"
        if length(couplings) != 3
            throw(ArgumentError("Incorrect number of couplings. Expected 2 (J,g), got $(length(couplings)). H = -J(XX + g1Z + g2X)"))
        # elseif couplings[2] == 0
        #     strings, coefficients = hamiltonian("TFIM", n, couplings[1:2], pbc=pbc)
        else
            strings = [generatexx(n, pbc=pbc) generatez(n) generatex(n)]
            coefficients = fill(-couplings[1], size(strings, 2))
            coefficients[end-2*n+1:end-n] *= couplings[2]
            coefficients[end-n+1:end] *= couplings[3]
        end

    elseif uppercase(model) == "HEISENBERG"
        length(couplings) == 1 || throw(ArgumentError("Incorrect number of couplings. Expected 1 (J), got $(length(couplings)). H = -J(XX + YY + ZZ)"))
        strings = [generatexx(n, pbc=pbc) generateyy(n, pbc=pbc) generatezz(n, pbc=pbc)]
        coefficients = fill(-couplings[1], size(strings, 2))

    elseif uppercase(model) == "XXZ"
        length(couplings) == 2 || throw(ArgumentError("Incorrect number of couplings. Expected 2 (J,Delta), got $(length(couplings)). H = -J(XX + YY + DeltaZZ)"))
        strings = [generatexx(n, pbc=pbc) generateyy(n, pbc=pbc) generatezz(n, pbc=pbc)]
        coefficients = fill(-couplings[1], size(strings, 2))
        coefficients[2*size(strings, 2)÷3+1:end] *= couplings[2]

elseif uppercase(model) == "GN"
        # length(couplings) == 2 || throw(ArgumentError("Incorrect number of couplings. Expected 2 (G,mu), got $(length(couplings)). H = (1+mu)(YX - XY) + GZ - GZZ "))
        nf > 0 || throw(ArgumentError("Number of flavors must be greater than 0."))
        strings = ones(Int8, n * nf, 0)
        for i in 0:nf-1
            temp = ones(Int8, n * nf, n - 1)
            temp[1+i*n:(i+1)*n, :] = generateyx(n)
            strings = [strings temp]
        end
        for i in 0:nf-1
            temp = ones(Int8, n * nf, n - 1)
            temp[1+i*n:(i+1)*n, :] = generatexy(n)
            strings = [strings temp]
        end
        strings = [strings generatez(n * nf)]
        for i in 0:nf-1
            for j in i+1:nf-1
                temp = ones(Int8, n * nf, n)
                temp[1+i*n:(i+1)*n, :] = generatez(n)
                temp[1+j*n:(j+1)*n, :] = generatez(n)
                strings = [strings temp]
            end
        end
        coefficients = [fill(1 - 2 * couplings[1], nf * (n - 1)); fill(-(1 - 2 * couplings[1]), nf * (n - 1)); fill(2 * couplings[3] - 4 * couplings[2] * (nf - 1), n * nf); fill(-4 * couplings[2], n * nf * (nf - 1) ÷ 2)]

    else
        throw(ArgumentError("Model not recognized."))
    end
    return strings, coefficients
end

function _dla(string1::AbstractVector{Int8}, iter)::Matrix{Int8}
    # result = Matrix{Int8}(undef, length(string1), 0)
    result = Vector{Int8}(undef, 0)
    for alg_string in iter
        string, _, c = pauliprod(string1, alg_string)
        # if !c
        #     string in eachcol(algebra) || (result = [result string])
        # end
        c || append!(result, string)
    end
    # return result
    return unique(reshape(result, length(string1), :), dims=2)
end

function dla(strings::AbstractMatrix{Int8})::Matrix{Int8}
    """
    dla
    ---
    This function generates the dynamical Lie algebra for a given set of Pauli strings.
    """
    n, finalind = size(strings)
    algebra = strings
    initialind = 1

    while true
        for i in axes(algebra, 2)[end:-1:initialind]
            iter = eachcol(algebra)[i-1:-1:1]
            chunks = Iterators.partition(iter, max(1, length(iter) ÷ Threads.nthreads()))
            tasks = map(chunks) do chunk
                Threads.@spawn _dla(view(algebra, :, i), chunk)
            end
            results = setdiff(fetch.(tasks), [Matrix{Int8}(undef, n, 0)])
            # for result in results
            #     result[:, 1] in algebra || (algebra = [algebra result[:, 1]])
            # end
            algebra = unique([algebra results...], dims=2)
        end
        finalind == size(algebra, 2) && return algebra
        initialind = finalind + 1
        finalind = size(algebra, 2)
    end
end

function subalgfind(strings::AbstractMatrix{Int8})::Matrix{Int8}
    """
    subalgfind
    ----------
    This function finds a Cartan subalgebra for a given set of Pauli strings.
    """
    sortedstrings = strings[:, sortperm(count(==(1), strings, dims=1)[1, :], rev=true)]
    subalgebra = sortedstrings[:, 1:1]
    for string1 in eachcol(sortedstrings)[2:end]
        for string2 in eachcol(subalgebra)
            if !(pauliprod(string1, string2)[3])
                break
            elseif string2 == @view subalgebra[:, end]
                subalgebra = [subalgebra string1]
                break
            end
        end
    end
    return subalgebra
end

function cartandecomp(strings::AbstractMatrix{Int8}, involution::Function)::Dict{String,Matrix{Int8}}
    """
    cartandecomp
    ------------
    This function finds the Cartan decomposition for a given set of Pauli strings.
    """
    inds = involution(strings)
    return Dict("k" => strings[:, inds], "m" => strings[:, .!inds], "h" => subalgfind(strings[:, .!inds]))
end