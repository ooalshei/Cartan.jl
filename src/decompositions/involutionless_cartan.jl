function _dla(string::Pair{<:Unsigned,<:Integer}, iter, Q::Integer)
    strings = Pair.(com.(string.first, first.(iter), Q), string.second .* last.(iter))
    return filter!(x -> !iszero(x.first), strings)
end

function involutionlessdecomp(paulis::PauliList)
    finalind = length(paulis)
    initialind = 1
    algebra = Pair.(paulis, Int8(-1))
    contradiction = false
    while !contradiction
        for i in eachindex(algebra)[end:-1:initialind]
            iter = @view algebra[i-1:-1:1]
            chunks = Iterators.partition(iter, max(1, length(iter) ÷ Threads.nthreads()))
            tasks = map(chunks) do chunk
                Threads.@spawn _dla(algebra[i], chunk, paulis.qubits)
            end
            union!(algebra, fetch.(tasks)...)
            length(algebra) == length(unique!(x -> x.first, algebra)) ||
                (contradiction = true)
        end
        if finalind == length(algebra)
            alg = PauliList(first.(algebra), paulis.qubits, iscopy=false)
            contradiction &&
                return Dict(:g => alg, :h => subalgfind(alg), :k => nothing, :m => nothing)
            m = filter(x -> x.second == -1, algebra)
            m = PauliList(first.(m), paulis.qubits, iscopy=false)
            k = filter(x -> x.second == 1, algebra)
            k = PauliList(first.(k), paulis.qubits, iscopy=false)
            return Dict(:g => alg, :k => k, :m => m, :h => subalgfind(m))
        end
        initialind = finalind + 1
        finalind = length(algebra)
        contradiction &&
            (algebra = PauliList(first.(algebra), paulis.qubits, iscopy=false); break)
    end

    while true
        for i in eachindex(algebra)[end:-1:initialind]
            iter = @view algebra[i-1:-1:1]
            chunks = Iterators.partition(iter, max(1, length(iter) ÷ Threads.nthreads()))
            tasks = map(chunks) do chunk
                Threads.@spawn _dla(algebra[i], chunk, paulis.qubits)
            end
            union!(algebra, fetch.(tasks)...)
        end

        finalind == length(algebra) && return Dict(
            :g => algebra,
            :h => subalgfind(algebra),
            :k => nothing,
            :m => nothing,
        )

        initialind = finalind + 1
        finalind = length(algebra)
    end
end
