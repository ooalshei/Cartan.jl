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
function subalgred(halg::PauliList)
    bstrings = PauliList(halg[1:1], halg.qubits, iscopy=false)
    multstrings = halg[1:1]
    i = 1
    while halg ⊈ multstrings
        for j in eachindex(halg)[i:end]
            i = j
            halg[j] ∈ multstrings || (push!(bstrings, halg[j]); break)
        end
        for string in multstrings
            push!(multstrings, string ⊻ halg[i])
        end
        push!(multstrings, halg[i])
    end
    return bstrings
end

"""
    fragmentedsubspaces(kstrings::PauliList, bstrings::PauliList)

Partition `kstrings` into subsets corresponding to each element of `bstrings` for a
reductive KHK decomposition.

For each `b` in `bstrings` the returned vector contains the
Pauli strings from `kstrings` that do not commute with `b` but commute with all the previous
elements in `bstrings`.

# Examples
```jldoctest
julia> ham = hamiltonian("TFIM", 4, [1, 1], UInt32)
PauliSentence{UInt32, Complex{Int64}, 4} with 7 entries:
  0x00000020 => -1+0im
  0x00000040 => -1+0im
  0x00000006 => -1+0im
  0x00000010 => -1+0im
  0x0000000c => -1+0im
  0x00000003 => -1+0im
  0x00000080 => -1+0im

julia> d = involutionlessdecomp(PauliList(collect(keys(ham)), 4))
Dict{Symbol, PauliList{UInt32, 4}} with 4 entries:
  :m => ["-Z--", "--Z-", "-XX-", "Z---", "--XX", "XX--", "---Z"…
  :k => ["--XY", "YX--", "XY--", "--YX", "-XY-", "-YX-", "-YZX"…
  :h => ["-Z--", "--Z-", "Z---", "---Z"]
  :g => ["-Z--", "--Z-", "-XX-", "Z---", "--XX", "XX--", "---Z"…

julia> println(d[:k])
["--XY", "YX--", "XY--", "--YX", "-XY-", "-YX-", "-YZX", "YZX-", "XZY-", "-XZY", "YZZX",
 "XZZY"]

julia> println(d[:h])
["-Z--", "--Z-", "Z---", "---Z"]

julia> fragmentedsubspaces(d[:k], d[:h])
4-element Vector{PauliList{UInt32, 4}}:
 ["YX--", "XY--", "-XY-", "-YX-", "-YZX", "-XZY"]
 ["--XY", "--YX", "YZX-", "XZY-"]
 ["YZZX", "XZZY"]
 0-element PauliList{UInt32, 4}
```
"""
function fragmentedsubspaces(
    kstrings::PauliList{T,Q},
    bstrings::PauliList{<:Unsigned,Q},
) where {T,Q}

    strings = copy(kstrings)
    symstrings = Vector{PauliList{T,Q}}(undef, length(bstrings))
    for i in eachindex(bstrings)
        temp = PauliList{T,Q}(undef, 0)
        j = 1
        while j <= length(strings)
            if com(strings[j], bstrings[i], Q).second
                j += 1
            else
                push!(temp, popat!(strings, j))
            end
        end
        symstrings[i] = temp
    end
    return symstrings
end

function cleangenerators!(symgenerators::AbstractVector{<:PauliList}, bstrings::PauliList)
    inds = findall(x -> length(x) == 0, symgenerators)
    deleteat!(symgenerators, inds)
    deleteat!(bstrings, inds)
end