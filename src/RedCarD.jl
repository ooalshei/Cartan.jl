"""
RedCarD

High-level utilities for Cartan decompositions and subalgebra
analysis on Symplectic Pauli algebras.

Purpose
- Provide builders and decomposition routines used for constructing and
    analyzing Cartan decompositions, involutions, and reductive
    subalgebras on Pauli-like operator algebras.

Exports
- `hamiltonian`, `algebra`, `dla`, `evenoddx`, `evenoddy`, `evenoddz`,
    `typeIorII`, `typeIII`, `subalgfind`, `cartandecomp`, `involutionlessdecomp`,
    `subalgred`, `fragmentedsubspaces`, `cleangenerators!`, `optimizer`,
    `reductive_optimizer`

Dependencies
- Re-exports selected symbols from `SymplecticPauli` (e.g. `Pauli`,
    `PauliList`, `PauliSentence`) and uses helper functions like `com`,
    `counti`, `countx`, `ad!`, etc.

Quick start
```
using RedCarD

# construct a Hamiltonian (example placeholder)
H = hamiltonian(PauliList([Pauli("X..")]))

# perform a Cartan decomposition
K, P = cartandecomp(H)
```

Place more detailed examples and API notes in the package README
or module-level docs as needed.
"""

module RedCarD

using Reexport
@reexport using SymplecticPauli: AbstractPauli, UPauli, Pauli, PauliList, PauliSentence
using SymplecticPauli: com, counti, countx, county, countz, ad, ad!

export hamiltonian,
    algebra,
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
