```@meta
CurrentModule = RedCarD
```

# API reference

```@docs
RedCarD
```

```@contents
Pages = ["api.md"]
Depth = 3
```

## Hamiltonian builders

```@docs
hamiltonian
```

## Algebra closure

```@docs
dla
subalgfind
```

## Involution maps

```@docs
evenoddx
evenoddy
evenoddz
typeIorII
typeIII
```

## Cartan splittings

```@docs
cartandecomp
involutionlessdecomp
```

## Reductive helpers

```@docs
subalgred
fragmentedsubspaces
cleangenerators!
```

## Angle solvers

```@docs
optimizer
reductive_optimizer
```

## Re-exported names

`using RedCarD` also brings in the Pauli representation these functions operate on:

| Name | Purpose |
|:--|:--|
| `UPauli`, `Pauli` | a single Pauli string, unsigned or with its ``\pm 1, \pm i`` phase |
| `PauliList` | an ordered list of strings — algebras, generator sets |
| `PauliSentence` | strings with coefficients — operators |
| `com` | commutation test, returning the product string and a commuting flag |
| `ad`, `ad!` | conjugation by ``\prod_j e^{i\theta_j g_j}`` |
| `tostring` | human-readable form (`"XY-Z"`) |
| `tomatrix` | dense matrix, for checking small cases |
| `countx`, `county`, `countz`, `counti` | bit counts on a string |

See [SymplecticPauli.jl](https://github.com/ooalshei/SymplecticPauli.jl) for their full
documentation.

## Index

```@index
```
