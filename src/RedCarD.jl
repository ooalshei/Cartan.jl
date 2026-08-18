@doc raw"""
    RedCarD

**Red**uctive **Car**tan **D**ecompositions of Pauli operator algebras, and the fixed-depth
circuits they produce.

A Hamiltonian ``H`` supported on a Lie algebra ``\mathfrak{g}`` that admits a Cartan
decomposition ``\mathfrak{g} = \mathfrak{k} \oplus \mathfrak{m}`` can be written as

```math
H = K h K^\dagger, \qquad K = \prod_j e^{i\theta_j k_j},\; k_j \in \mathfrak{k},
\quad h \in \mathfrak{h},
```

with ``\mathfrak{h} \subseteq \mathfrak{m}`` a Cartan subalgebra. Because everything in
``\mathfrak{h}`` commutes, ``e^{-iHt} = K e^{-iht} K^\dagger`` is an *exact* circuit whose
depth does not grow with ``t`` — the whole point of the construction.

This package covers that pipeline end to end:

| Step | Functions |
|:--|:--|
| Build a model Hamiltonian | [`hamiltonian`](@ref) |
| Close the algebra it generates | [`dla`](@ref) |
| Split it | [`cartandecomp`](@ref), [`involutionlessdecomp`](@ref) |
| Find a Cartan subalgebra | [`subalgfind`](@ref), [`subalgred`](@ref), [`fragmentedsubspaces`](@ref) |
| Solve for the angles | [`optimizer`](@ref), [`reductive_optimizer`](@ref) |

Pauli operators, their symplectic bit-string representation and the rotations acting on them
come from [`SymplecticPauli`](https://github.com/ooalshei/SymplecticPauli.jl), whose core
API (`PauliList`, `PauliSentence`, `com`, `ad`, `tostring`, `tomatrix`, …) is re-exported, so
`using RedCarD` is enough.

# Examples
```jldoctest
julia> H = hamiltonian("TFIM", 4, [1.0, 0.5]);     # -J Σ XᵢXᵢ₊₁ - Jg Σ Zᵢ

julia> decomposition = involutionlessdecomp(PauliList(collect(keys(H)), 4));

julia> length(decomposition.g), length(decomposition.k), length(decomposition.m)
(28, 12, 16)

julia> tostring(decomposition.h)     # the Cartan subalgebra to rotate into
4-element Vector{String}:
 "-Z--"
 "--Z-"
 "Z---"
 "---Z"
```

The algorithms follow *Cartan decompositions for Pauli operator algebras*,
[arXiv:2512.06070](https://arxiv.org/abs/2512.06070).
"""
module RedCarD

using Reexport
@reexport using SymplecticPauli:
    AbstractPauli,
    UPauli,
    Pauli,
    PauliList,
    PauliSentence,
    ad,
    ad!,
    com,
    counti,
    countx,
    county,
    countz,
    tomatrix,
    tostring
using SymplecticPauli: _check_string_length, _ipow, _prodsign, _symplectic_prod

export hamiltonian,
    dla,
    evenoddx,
    evenoddy,
    evenoddz,
    typeIorII,
    typeIII,
    subalgfind,
    cartandecomp,
    involutionlessdecomp,
    subalgred,
    fragmentedsubspaces,
    cleangenerators!,
    optimizer,
    reductive_optimizer

include("builders.jl")
include("decompositions/involutions.jl")
include("decompositions/cartan.jl")
include("decompositions/involutionless_cartan.jl")
include("decompositions/reductive_cartan.jl")
include("optimizers/optimizer.jl")
include("optimizers/reductive_optimizer.jl")

end
