"""
    subalgred(halg::PauliList)

Return a minimal set of independent Pauli strings whose multiplicative closure (under `⊻`)
generates the algebra spanned by `halg`.

Perform a greedy construction: start from the first string in `halg` and iteratively pick
the next element not already in the set of generated products, adding it as a generator and
updating the generated products by multiplying the new generator with the existing products.

# Examples
```jldoctest
julia> h = PauliList(["Z--", "-Z-", "ZZ-", "--Z", "Z-Z", "-ZZ", "ZZZ"])
7-element PauliList{UInt64, 3}:
 0x0000000000000008
 0x0000000000000010
 0x0000000000000018
 0x0000000000000020
 0x0000000000000028
 0x0000000000000030
 0x0000000000000038

julia> b = subalgred(h)
3-element PauliList{UInt64, 3}:
 0x0000000000000008
 0x0000000000000010
 0x0000000000000020

julia> print(b)
["Z--", "-Z-", "--Z"]
```
"""
function subalgred(halg::PauliList{T,Q}) where {T,Q}
    isempty(halg) && return PauliList{T,Q}(T[], iscopy=false, check=false)
    bstrings = T[halg[1]]
    # Products generated so far, held in a set: the closure doubles in size with every new
    # generator, so testing membership against a list would dominate the construction.
    multstrings = Set{T}((halg[1],))
    i = 1
    while !all(in(multstrings), halg)
        for j in eachindex(halg)[i:end]
            i = j
            halg[j] ∈ multstrings || (push!(bstrings, halg[j]); break)
        end
        for string in collect(multstrings)
            push!(multstrings, string ⊻ halg[i])
        end
        push!(multstrings, halg[i])
    end
    return PauliList{T,Q}(bstrings, iscopy=false, check=false)
end

"""
    fragmentedsubspaces(kstrings::PauliList, bstrings::PauliList)

Partition `kstrings` into subsets corresponding to each element of `bstrings` for a
reductive KHK decomposition.

For each `b` in `bstrings` the returned vector contains the
Pauli strings from `kstrings` that do not commute with `b` but commute with all the previous
elements in `bstrings`.

See also [`subalgred`](@ref), [`cleangenerators!`](@ref), [`reductive_optimizer`](@ref).

# Examples
```jldoctest
julia> ham = hamiltonian("TFIM", 4, [1, 1], UInt32);

julia> decomposition = involutionlessdecomp(PauliList(collect(keys(ham)), 4));

julia> tostring(decomposition.h)
4-element Vector{String}:
 "-Z--"
 "--Z-"
 "Z---"
 "---Z"

julia> fragments = fragmentedsubspaces(decomposition.k, decomposition.h);

julia> length.(fragments)
4-element Vector{Int64}:
 6
 4
 2
 0

julia> tostring(fragments[1])       # the k elements that move the first Z
6-element Vector{String}:
 "YX--"
 "XY--"
 "-XY-"
 "-YX-"
 "-YZX"
 "-XZY"
```
"""
function fragmentedsubspaces(
    kstrings::PauliList{T,Q},
    bstrings::PauliList{<:Unsigned,Q},
) where {T,Q}

    strings = copy(kstrings.strings)
    symstrings = Vector{PauliList{T,Q}}(undef, length(bstrings))
    for i in eachindex(bstrings)
        # Splitting into two fresh vectors keeps this linear: deleting the matches from
        # `strings` one at a time shifts the tail on every hit.
        anticommuting = T[]
        remaining = T[]
        for string in strings
            if com(string, bstrings[i], Q).second
                push!(remaining, string)
            else
                push!(anticommuting, string)
            end
        end
        symstrings[i] = PauliList{T,Q}(anticommuting, iscopy=false, check=false)
        strings = remaining
    end
    return symstrings
end

"""
    cleangenerators!(symgenerators::AbstractVector{<:PauliList}, bstrings::PauliList)

Drop the empty fragments from `symgenerators`, together with the elements of `bstrings` they
belong to.

A Cartan generator with no ``\\mathfrak{k}`` elements acting on it contributes nothing to the
optimization, and [`reductive_optimizer`](@ref) expects the two arguments to line up, so
both are shortened in place.

See also [`fragmentedsubspaces`](@ref).

# Examples
```jldoctest
julia> ham = hamiltonian("TFIM", 4, [1, 1], UInt32);

julia> decomposition = involutionlessdecomp(PauliList(collect(keys(ham)), 4));

julia> bstrings = subalgred(decomposition.h);

julia> fragments = fragmentedsubspaces(decomposition.k, bstrings);

julia> length.(fragments)
4-element Vector{Int64}:
 6
 4
 2
 0

julia> cleangenerators!(fragments, bstrings);

julia> length.(fragments), length(bstrings)
([6, 4, 2], 3)
```
"""
function cleangenerators!(symgenerators::AbstractVector{<:PauliList}, bstrings::PauliList)
    inds = findall(x -> length(x) == 0, symgenerators)
    deleteat!(symgenerators, inds)
    deleteat!(bstrings, inds)
end
