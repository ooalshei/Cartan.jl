_dla(string::Unsigned, iter, Q::Integer) =
    PauliList(filter!(!iszero, com.(string, iter, Q)), Q, iscopy=false)

function dla(paulis::PauliList)
    finalind = length(paulis)
    initialind = 1
    algebra = copy(paulis)
    while true
        for i in eachindex(algebra)[end:-1:initialind]
            iter = @view algebra[i-1:-1:1]
            chunks = Iterators.partition(iter, max(1, length(iter) ÷ Threads.nthreads()))
            tasks = map(chunks) do chunk
                Threads.@spawn _dla(algebra[i], chunk, paulis.qubits)
            end
            union!(algebra, fetch.(tasks)...)
        end
        finalind == length(algebra) && return algebra
        initialind = finalind + 1
        finalind = length(algebra)
    end
end

function subalgfind(paulis::PauliList)
    sortedstrings = sort(paulis, by=x -> counti(x, paulis.qubits), rev=true)
    subalgebra = sortedstrings[1:1]
    for g in sortedstrings[2:end]
        for h in subalgebra
            iszero(com(g, h, paulis.qubits)) || break
            h == subalgebra[end] && (push!(subalgebra, g); break)
        end
    end
    return PauliList(subalgebra, paulis.qubits, iscopy=false)
end

function cartandecomp(paulis::PauliList, involution::Function)
    inds = involution(paulis)
    k = PauliList(paulis[inds], paulis.qubits, iscopy=false)
    m = PauliList(paulis[.!inds], paulis.qubits, iscopy=false)
    return Dict(:k => k, :m => m, :h => subalgfind(m))
end
