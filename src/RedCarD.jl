"""
    RedCarD.jl

High-level documentation for the RedCarD.jl package.

This package implements the Cartan decomposition algorithm for Pauli operator
algebras described in "Cartan decompositions for Pauli operator algebras"
(arXiv:2512.06070). The algorithm used by the package is summarized on page 5
of the paper; the repository implements the routines to build Pauli-based
Hamiltonians, compute decompositions (including involutionless and reductive
variants), and perform the roto-optimization to block-diagonalize Hamiltonians
with respect to an abelian subalgebra.

Usage
    - Use `hamiltonian` to build common lattice models (Ising, XY, TFIM, ...).
    - Use `dla`, `cartandecomp`, or `involutionlessdecomp` to compute decompositions.
    - Use `optimizer` or `reductive_optimizer` to rotate a Hamiltonian into a
        chosen subalgebra according to the algorithm.

See also the package README for examples and the source paper for algorithmic
details: https://arxiv.org/pdf/2512.06070
"""
module RedCarD

using Reexport
@reexport using SymplecticPauli: AbstractPauli, UPauli, Pauli, PauliList, PauliSentence
using SymplecticPauli: com, counti, countx, county, countz, ad, ad!, _check_string_length

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
