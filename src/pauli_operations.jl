"""
Pauli Operations
----------------
This module contains several useful operations on Pauli strings. The identity along with the three Pauli matrices are
interchangeably referred to as (1, 2, 3, 4) or (-, X, Y, Z).
"""

# These arrays are used to find products of Pauli matrices.
const rules = Int64[1 2 3 4
    2 1 4 3
    3 4 1 2
    4 3 2 1]
const sign_rules = ComplexF64[1 1 1 1
    1 1 1im -1im
    1 -1im 1 1im
    1 1im -1im 1]
print(sign_rules[4, 3])


function product(string1::Vector{Int64}, string2::Vector{Int64})
    if length(string1) != length(string2)
        throw("Dimension mismatch ($(length(string1)) and $(length(string2)))")
    end

    index = CartesianIndex.(string1, string2)
    sign = prod(sign_rules[index])
    result = rules[index]
    if imag(sign) == 0
        return result, sign, true
    else
        return result, sign, false
    end
end

x = [1, 2, 3]
y = [2, 3, 3]
product(x, y)
