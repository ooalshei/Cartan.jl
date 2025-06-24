module Cartan

export pauliprod, paulisum, paulidiff, pauliprod, conjugate, hamiltonian, algebra, dla, evenoddx, evenoddy, evenoddz, subalgfind, cartandecomp, involutionlessdecomp, subalgred, symsubspaces, cleangenerators!, optimizer, iterativeoptimizer, Symplectic

include("pauli_operations.jl")
include("decompositions/involutions.jl")
include("decompositions/cartan.jl")
include("decompositions/involutionless_cartan.jl")
include("decompositions/iterative_cartan.jl")
include("optimizers/optimizer.jl")
include("optimizers/iterative_optimizer.jl")

module Symplectic
# https://github.com/ooalshei/SymplecticPauli.jl/tree/dev
include("../../SymplecticPauli.jl/src/SymplecticPauli.jl")
using ..Cartan: _mutirr, _minanglefind
using .SymplecticPauli

include("symplectic_optimizers/optimizer.jl")
include("symplectic_optimizers/iterative_optimizer.jl")

end

end
