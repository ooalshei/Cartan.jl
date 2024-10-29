"""
cartan
------
This module contains methods that build Pauli strings given a Hamiltonian, generate the dynamical Lie algebra, and find
Cartan decompositions. The identity along with the three Pauli matrices are referred to as (1, 2, 3, 4).
"""
function generatexx(n::Int; pbc::Bool=false)::Matrix{Int8}
    """
    generatexx
    ----------
    This function generates the XX strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] .= 2
    end
    if pbc
        strings = [strings ones(Int8, n)]
        strings[[n, 1], n] .= 2
    end
    return strings
end

function generateyy(n::Int; pbc::Bool=false)::Matrix{Int8}
    """
    generateyy
    ----------
    This function generates the YY strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] .= 3
    end
    if pbc
        strings = [strings ones(Int8, n)]
        strings[[n, 1], n] .= 3
    end
    return strings
end

function generatezz(n::Int; pbc::Bool=false)::Matrix{Int8}
    """
    generatezz
    ----------
    This function generates the ZZ strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] .= 4
    end
    if pbc
        strings = [strings ones(Int8, n)]
        strings[[n, 1], n] .= 4
    end
    return strings
end

function generatexy(n::Int)::Matrix{Int8}
    """
    generatexy
    ----------
    This function generates the XY strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] = [2, 3]
    end
    return strings
end

function generateyx(n::Int)::Matrix{Int8}
    """
    generateyx
    ----------
    This function generates the YX strings for a chain of length n.
    """
    strings = ones(Int8, n, n - 1)
    for i in 1:n-1
        strings[[i, i + 1], i] = [3, 2]
    end
    return strings
end

function generatex(n::Int)::Matrix{Int8}
    """
    generatex
    ---------
    This function generates the X strings for a chain of length n.
    """
    strings = ones(Int8, n, n)
    strings[CartesianIndex.(1:n, 1:n)] .= 2
    return strings
end

function generatey(n::Int)::Matrix{Int8}
    """
    generatey
    ---------
    This function generates the Y strings for a chain of length n.
    """
    strings = ones(Int8, n, n)
    strings[CartesianIndex.(1:n, 1:n)] .= 3
    return strings
end

function generatez(n::Int)::Matrix{Int8}
    """
    generatez
    ---------
    This function generates the Z strings for a chain of length n.
    """
    strings = ones(Int8, n, n)
    strings[CartesianIndex.(1:n, 1:n)] .= 4
    return strings
end

function hamiltonian(model::String, n::Int, couplings::Vector{Float64}=[1.0]; pbc::Bool=false)::Tuple{Matrix{Int8},Vector{Float64}}
    """
    hamiltonian
    -----------
    This function generates the Hamiltonian for a given model.
    """
    if uppercase(model) == "ISING"
        length(couplings) == 1 || ArgumentError("Incorrect number of couplings. Expected 1 (J), got $(length(couplings)). H = -JXX")
        strings = generatexx(n, pbc=pbc)
        coefficients = -couplings[1] * ones(Float64, size(strings, 2))

    elseif uppercase(model) == "XY"
        length(couplings) == 1 || ArgumentError("Incorrect number of couplings. Expected 1 (J), got $(length(couplings)). H = -J(XX + YY)")
        strings = [generatexx(n, pbc=pbc) generateyy(n, pbc=pbc)]
        coefficients = -couplings[1] * ones(Float64, size(strings, 2))

    elseif uppercase(model) == "TFIM"
        if length(couplings) != 2
            ArgumentError("Incorrect number of couplings. Expected 2 (J,g), got $(length(couplings)). H = -J(XX + gZ)")
        elseif couplings[2] == 0
            strings, coefficients = hamiltonian("Ising", n, [couplings[1]], pbc=pbc)
        else
            strings = [generatexx(n, pbc=pbc) generatez(n)]
            coefficients = -couplings[1] * ones(Float64, size(strings, 2))
            coefficients[end-n+1:end] *= couplings[2]
        end

    elseif uppercase(model) == "TFXY"
        if length(couplings) != 2
            ArgumentError("Incorrect number of couplings. Expected 2 (J,g), got $(length(couplings)). H = -J(XX + YY + gZ)")
        elseif couplings[2] == 0
            strings, coefficients = hamiltonian("XY", n, [couplings[1]], pbc=pbc)
        else
            strings = [generatexx(n, pbc=pbc) generateyy(n, pbc=pbc) generatez(n)]
            coefficients = -couplings[1] * ones(Float64, size(strings, 2))
            coefficients[end-n+1:end] *= couplings[2]
        end

    elseif uppercase(model) == "HEISENBERG"
        length(couplings) == 1 || ArgumentError("Incorrect number of couplings. Expected 1 (J), got $(length(couplings)). H = -J(XX + YY + ZZ)")
        strings = [generatexx(n, pbc=pbc) generateyy(n, pbc=pbc) generatezz(n, pbc=pbc)]
        coefficients = -couplings[1] * ones(Float64, size(strings, 2))

    elseif uppercase(model) == "GN"
        length(couplings) == 1 || ArgumentError("Incorrect number of couplings. Expected 2 (G,mu), got $(length(couplings)). H = (1+mu)(YX - XY) + GZ - GZZ ")
        strings = [generateyx(n) generatexy(n) generatez(n) [generatez(div(n, 2)); generatez(div(n, 2))]]
        coefficients = [(0.5 + couplings[1]) * ones(Float64, n - 1); -(0.5 + couplings[1]) * ones(Float64, n - 1); couplings[2] * ones(Float64, n); -couplings[2] * ones(Float64, div(n, 2))]

    else
        ArgumentError("Model not recognized.")
    end
    return strings, coefficients
end

function _dla(algebra::Matrix{Int8}, string1::Vector{Int8}, iter)::Matrix{Int8}
    result = Matrix{Int8}(undef, length(string1), 0)
    for j in iter
        string, _, c = pauliprod(string1, algebra[:, j])
        if !c
            string in eachcol(algebra) || (result = [result string])
        end
    end
    return result
end

function dla(strings::Matrix{Int8})::Matrix{Int8}
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
            iter = axes(algebra, 2)[i-1:-1:1]
            chunks = Iterators.partition(iter, max(1, length(iter) ÷ Threads.nthreads()))
            tasks = map(chunks) do chunk
                Threads.@spawn _dla(algebra, algebra[:, i], chunk)
            end
            results = setdiff(fetch.(tasks), [Matrix{Int8}(undef, n, 0)])
            for result in results
                result[:, 1] in algebra || (algebra = [algebra result[:, 1]])
            end
        end
        finalind == size(algebra, 2) && return algebra
        initialind = finalind + 1
        finalind = size(algebra, 2)
    end
end

function subalgfind(strings::Matrix{Int8})::Matrix{Int8}
    """
    subalgfind
    ----------
    This function finds a Cartan subalgebra for a given set of Pauli strings.
    """
    sortedstrings = strings[:, sortperm(count(==(1), strings, dims=1)[1, :], rev=true)]
    subalgebra = sortedstrings[:, 1:1]
    for i in axes(sortedstrings, 2)[2:end]
        for j in axes(subalgebra, 2)
            if !(pauliprod(sortedstrings[:, i], subalgebra[:, j])[3])
                break
            elseif j == size(subalgebra, 2)
                subalgebra = [subalgebra sortedstrings[:, i]]
                break
            end
        end
    end
    return subalgebra
end

function cartandecomp(strings::Matrix{Int8}, involution::Function)::Dict{String,Matrix{Int8}}
    """
    cartandecomp
    ------------
    This function finds the Cartan decomposition for a given set of Pauli strings.
    """
    inds = involution(strings)
    return Dict("k" => strings[:, inds], "m" => strings[:, .!inds], "h" => subalgfind(strings[:, .!inds]))
end

