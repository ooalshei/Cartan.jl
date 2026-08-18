```@meta
CurrentModule = RedCarD
```

# Algebras and decompositions

## Pauli strings

Everything here is a symplectic bit string: `2Q` bits per operator on `Q` qubits, the low
half recording `X` support and the high half `Z` support, so a `Y` sets both bits of its
qubit. Products are `⊻`, and whether two strings commute is one popcount.

Two containers carry them:

* `PauliList` — an ordered list of strings, used for algebras and generator sets.
* `PauliSentence` — strings with coefficients, used for operators.

```jldoctest algebras
julia> using RedCarD

julia> paulis = PauliList(["XY-", "-ZZ"]);   # stored as bit strings

julia> tostring(paulis)                     # …and read back as text
2-element Vector{String}:
 "XY-"
 "-ZZ"

julia> com(paulis[1], paulis[2], 3).second      # do they commute?
false

julia> H = PauliSentence(paulis, [0.5, -1.0]);

julia> length(H)
2
```

## Building Hamiltonians

[`hamiltonian`](@ref) covers the usual spin chains, with open or periodic boundaries.

```jldoctest algebras
julia> length(hamiltonian("HEISENBERG", 4, [1.0]))         # 3 bonds × 3 terms
9

julia> length(hamiltonian("HEISENBERG", 4, [1.0], pbc=true))
12

julia> sort(collect(keys(tostring(hamiltonian("XXZ", 3, [1.0, 0.4])))))
6-element Vector{String}:
 "-XX"
 "-YY"
 "-ZZ"
 "XX-"
 "YY-"
 "ZZ-"
```

The narrower the integer type, the less memory and the faster the bit work; `UInt32` holds
up to 16 qubits and is worth passing explicitly for large runs:

```jldoctest algebras
julia> H32 = hamiltonian("TFIM", 12, [1.0, 0.5], UInt32);

julia> keytype(H32)
UInt32
```

## Closing the algebra

[`dla`](@ref) takes commutators until nothing new appears.

```jldoctest algebras
julia> closure = dla(PauliList(["X-", "-X", "Z-", "-Z"]));

julia> tostring(closure)
6-element Vector{String}:
 "X-"
 "-X"
 "Z-"
 "-Z"
 "-Y"
 "Y-"
```

The dimension of the closure decides whether any of this is worth doing: an algebra that
grows to $4^Q$ gives a circuit no shallower than the generic one. Check it before
committing to a decomposition — closures are cheap compared to what follows.

## Decomposition by involution

A Cartan decomposition can be produced from an involution $\theta$ — a map with
$\theta^2 = 1$ that respects the bracket — by taking $\mathfrak{k}$ and $\mathfrak{m}$ to be
its $+1$ and $-1$ eigenspaces. [`cartandecomp`](@ref) takes the involution as a function
returning one boolean per element.

Three families are provided:

| Involution | Rule | Functions |
|:--|:--|:--|
| Even-odd counting | parity of the number of `X`/`Y`/`Z` factors | [`evenoddx`](@ref), [`evenoddy`](@ref), [`evenoddz`](@ref) |
| Type I/II | $-q p^{\mathsf T} q = \pm p$ for a fixed string `q` | [`typeIorII`](@ref) |
| Type III | $q p q = \pm p$ for a fixed string `q` | [`typeIII`](@ref) |

```jldoctest algebras
julia> split = cartandecomp(closure, evenoddy);

julia> tostring(split.k), tostring(split.m)
(["-Y", "Y-"], ["X-", "-X", "Z-", "-Z"])

julia> tostring(split.h)
2-element Vector{String}:
 "X-"
 "-X"
```

The conjugating-string involutions need the string supplied, so wrap them:

```jldoctest algebras
julia> other = cartandecomp(closure, p -> typeIII(p, UPauli("ZZ")));

julia> tostring(other.k)
2-element Vector{String}:
 "Z-"
 "-Z"
```

A decomposition is only useful when the Hamiltonian lands entirely in $\mathfrak{m}$, which
is a property of the involution *and* of the model. Picking one that works is the reason for
the next section.

## Decomposition without an involution

[`involutionlessdecomp`](@ref) sidesteps the search. It closes the algebra while carrying a
sign per element — $-1$ for everything the Hamiltonian generates, and the product of signs
through every commutator — which is exactly the condition
$[\mathfrak{m},\mathfrak{m}] \subseteq \mathfrak{k}$ written out. If no string is ever
reached with both signs, the result is a genuine Cartan decomposition that keeps the
Hamiltonian in $\mathfrak{m}$ by construction.

```jldoctest algebras
julia> H = hamiltonian("TFIM", 4, [1.0, 0.5]);

julia> decomposition = involutionlessdecomp(PauliList(collect(keys(H)), 4));

julia> length(decomposition.k), length(decomposition.m)
(12, 16)
```

When a contradiction *is* found, no such decomposition exists; `k` and `m` come back as
`nothing`, and `g` and `h` still describe the closed algebra and a commuting subalgebra of
it:

```jldoctest algebras
julia> failed = involutionlessdecomp(PauliList(["X-", "Y-", "Z-"]));

julia> failed.k === nothing, failed.m === nothing
(true, true)

julia> length(failed.g)
3
```

## Reductive decompositions

The reductive scheme rotates into the Cartan subalgebra one generator at a time rather than
all at once. It needs two things: a multiplicative generating set for $\mathfrak{h}$, and
the pieces of $\mathfrak{k}$ that act on each of them.

[`subalgred`](@ref) reduces $\mathfrak{h}$ to the strings whose products span it:

```jldoctest algebras
julia> h = PauliList(["Z--", "-Z-", "ZZ-", "--Z", "Z-Z", "-ZZ", "ZZZ"]);

julia> tostring(subalgred(h))
3-element Vector{String}:
 "Z--"
 "-Z-"
 "--Z"
```

[`fragmentedsubspaces`](@ref) then splits $\mathfrak{k}$ by which generator each element
first fails to commute with, and [`cleangenerators!`](@ref) drops the generators that end up
with nothing to act on:

```jldoctest algebras
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

These are exactly the inputs [`reductive_optimizer`](@ref) expects.
