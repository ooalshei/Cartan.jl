"""
This file contains different involutions.
"""
# evenoddx(strings::PauliList) = iseven.(countx.(strings))
# evenoddy(strings::PauliList) = iseven.(county.(strings))
# evenoddz(strings::PauliList) = iseven.(countz.(strings))

for p in [:x, :y, :z]
    @eval $(Symbol(:evenodd, p))(paulis::PauliList) = isodd.(count$(p).(paulis))
end

function typeIorII(paulis::PauliList, string::Unsigned)
    _check_string_length(string, paulis.qubits)
    result = Vector{Bool}(undef, length(paulis))
    for (i, pauli) in enumerate(paulis)
        isequal(com(pauli, string, paulis.qubits), 0) ? (result[i] = isodd(county(pauli))) : (result[i] = iseven(county(pauli)))
    end
    return result
end
typeIorII(paulis::PauliList{<:Unsigned,Q}, string::AbstractPauli{<:Unsigned,Q}) where {Q} = typeIorII(paulis, string.string)

function typeIII(paulis::PauliList, string::Unsigned)
    result = Vector{Bool}(undef, length(paulis))
    for (i, pauli) in enumerate(paulis)
        isequal(com(pauli, string, paulis.qubits), 0) ? (result[i] = 1) : (result[i] = 0)
    end
    return result
end
typeIII(paulis::PauliList{<:Unsigned,Q}, string::AbstractPauli{<:Unsigned,Q}) where {Q} = typeIII(paulis, string.string)