function subalgred(halg::PauliList)
    abstrings = PauliList(halg[1:1], halg.qubits, iscopy=false)
    multstrings = halg[1:1]
    i = 1
    while halg ⊈ multstrings
        for j in eachindex(halg)[i:end]
            i = j
            halg[j] ∈ multstrings || (push!(abstrings, halg[j]); break)
        end
        for string in multstrings
            push!(multstrings, string ⊻ halg[i])
        end
        push!(multstrings, halg[i])
    end
    return abstrings
end

function symsubspaces(kstrings::PauliList{T,Q}, abstrings::PauliList{<:Unsigned,Q}) where {T,Q}
    strings = copy(kstrings)
    symstrings = Vector{PauliList{T,Q}}(undef, length(abstrings))
    for i in eachindex(abstrings)
        temp = PauliList{T,Q}(undef, 0)
        j = 1
        while j <= length(strings)
            iszero(com(strings[j], abstrings[i], Q)) ? j += 1 : push!(temp, popat!(strings, j))
        end
        symstrings[i] = temp
    end
    return symstrings
end

function cleangenerators!(symgenerators::AbstractVector{<:PauliList}, abstrings::PauliList)

    inds = findall(x -> length(x) == 0, symgenerators)
    deleteat!(symgenerators, inds)
    deleteat!(abstrings, inds)
end