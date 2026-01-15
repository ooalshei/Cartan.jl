for p in [:x, :y, :z]
    @eval begin
        """
            $($(Symbol(:evenodd, p)))

        Even-odd parity of the number of `$($(uppercase(string(p))))`'s in each string.
        """
        $(Symbol(:evenodd, p))(paulis::PauliList) =
            isodd.($(Symbol(:count, p)).(paulis, paulis.qubits))
    end
end

@doc raw"""
    typeIorII(p, q)

Type I or II involution. Find if ``-qp^{T}q = ± p``.
"""
function typeIorII end

function typeIorII(paulis::PauliList, string::Unsigned)
    _check_string_length(string, paulis.qubits)
    result = Vector{Bool}(undef, length(paulis))
    for (i, pauli) in enumerate(paulis)
        if com(pauli, string, paulis.qubits).second
            (result[i] = isodd(county(pauli, paulis.qubits)))
        else
            (result[i] = iseven(county(pauli, paulis.qubits)))
        end
    end
    return result
end
typeIorII(paulis::PauliList{<:Unsigned,Q}, string::AbstractPauli{<:Unsigned,Q}) where {Q} =
    typeIorII(paulis, string.string)

@doc raw"""
    typeIII(p, q)

Type III involution. Find if ``qpq = ± p``.
"""
function typeIII end

function typeIII(paulis::PauliList, string::Unsigned)
    result = Vector{Bool}(undef, length(paulis))
    for (i, pauli) in enumerate(paulis)
        com(pauli, string, paulis.qubits).second ? (result[i] = 1) : (result[i] = 0)
    end
    return result
end
typeIII(paulis::PauliList{<:Unsigned,Q}, string::AbstractPauli{<:Unsigned,Q}) where {Q} =
    typeIII(paulis, string.string)