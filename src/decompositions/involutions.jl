"""
This file contains different involutions.
"""
function evenoddx(strings::AbstractMatrix{Int8})::Vector{Bool}
    """
    evenoddx
    --------
    This function finds the even-odd X Pauli count for a given set of Pauli strings.
    """
    return ((count(==(2), strings, dims=1) .% 2))[1, :]
end

function evenoddy(strings::AbstractMatrix{Int8})::Vector{Bool}
    """
    evenoddy
    --------
    This function finds the even-odd Y Pauli count for a given set of Pauli strings.
    """
    return (count(==(3), strings, dims=1) .% 2)[1, :]
end

function evenoddz(strings::AbstractMatrix{Int8})::Vector{Bool}
    """
    evenoddz
    --------
    This function finds the even-odd Z Pauli count for a given set of Pauli strings.
    """
    return (count(==(4), strings, dims=1) .% 2)[1, :]
end