module Cartan

@static if Sys.isunix()
    run(`bash .install_git.sh`)
elseif Sys.iswindows()
    run(`cmd .install_git.cmd`)  
end
run(`git submodule update --init`)

using Reexport
include("../SymplecticPauli.jl/src/SymplecticPauli.jl")
@reexport using .SymplecticPauli: AbstractPauli, UPauli, Pauli, PauliList, PauliSentence
using .SymplecticPauli: com, counti, countx, county, countz, ad, ad!

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
    symsubspaces,
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
