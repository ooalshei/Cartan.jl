raw"""
involutionless_cartan
------
This module finds a Cartan decomposition such that :math:`\mathcal{H} \subset \mathfrak{m}` without the need for an
involution. This decomposition is unique, and if no such decomposition is found, the Hamiltonian is said to be
non-Cartan. The identity along with the three Pauli matrices are referred to as (1, 2, 3, 4).
"""
const flag_rules = Int8[1 2 3 4 0 0;
    2 1 4 3 0 0;
    3 4 1 2 0 0;
    4 3 2 1 0 0;
    0 0 0 0 6 5;
    0 0 0 0 5 6]

const flag_signs = Complex{Int8}[1 1 1 1 0 0;
    1 1 1im -1im 0 0;
    1 -1im 1 1im 0 0;
    1 1im -1im 1 0 0;
    0 0 0 0 1 1;
    0 0 0 0 1 1]

function _flag_pauliprod(string1::Vector{Int8},
    string2::Vector{Int8})::Tuple{Vector{Int8},Complex{Int8},Bool}

    length(string1) == length(string2) || throw(DimensionMismatch("Strings have different lengths ($(length(string1)) and $(length(string2)))"))
    index = CartesianIndex.(string1, string2)
    result = flag_rules[index]
    sign = prod(flag_signs[index])
    imag(sign) == 0 ? (return result, sign, true) : (return result, sign, false)
end

function _flag_dla(algebra::Matrix{Int8}, string1::Vector{Int8}, iter)::Matrix{Int8}
    # result = Matrix{Int8}(undef, length(string1), 0)
    result = Vector{Int8}(undef, 0)
    for j in iter
        string, _, c = _flag_pauliprod(string1, algebra[:, j])
        # c || (result = [result string])
        c || append!(result, string)
    end
    return unique(reshape(result, length(string1), :), dims=2)
end

function involutionlessdecomp(strings::Matrix{Int8})::Dict
    """
    dla
    ---
    This function generates the dynamical Lie algebra for a given set of Pauli strings. It performs the Cartan decomposition that places the generating set in m, if it exists. Otherwise, it finds a Cartan subalgebra within the entire dla.
    """
    n, finalind = size(strings)
    algebra = Int8[strings; 5 * ones(Int8, 1, size(strings, 2))]
    initialind = 1
    contradiction = false
    while true
        for i in axes(algebra, 2)[end:-1:initialind]
            j = axes(algebra, 2)[i-1:-1:1]
            chunks = Iterators.partition(j, max(1, length(j) ÷ Threads.nthreads()))
            tasks = map(chunks) do chunk
                Threads.@spawn _flag_dla(algebra, algebra[:, i], chunk)
            end
            tasks = fetch.(tasks)
            flag_results = setdiff(tasks, [Matrix{Int8}(undef, n + 1, 0)])
            flag_results = unique(hcat(flag_results...), dims=2)
            results = unique(flag_results[1:end-1, :], dims=2)
            size(results, 2) == size(flag_results, 2) || (contradiction = true)
            size(flag_results, 1) == size(algebra, 1) && (algebra = unique([algebra flag_results], dims=2))
            
        end

        if finalind == size(algebra, 2)
            if contradiction
                return Dict("g" => algebra[1:end-1, :], "h" => subalgfind(algebra[1:end-1, :]), "k" => nothing, "m" => nothing)
            else
                ind5 = findall(==(5), algebra[end, :])
                ind6 = findall(==(6), algebra[end, :])
                mset = algebra[1:end-1, ind5]
                return Dict("g" => algebra,"k" => algebra[1:end-1, ind6], "m" => mset, "h" => subalgfind(mset))
            end
        end
        initialind = finalind + 1
        finalind = size(algebra, 2)
    end

end
