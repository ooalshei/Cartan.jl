# RedCarD.jl

[![Build Status](https://github.com/ooalshei/RedCarD.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ooalshei/RedCarD.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ooalshei/RedCarD.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ooalshei/RedCarD.jl)
[![Coverage](https://coveralls.io/repos/github/ooalshei/RedCarD.jl/badge.svg?branch=main)](https://coveralls.io/github/ooalshei/RedCarD.jl?branch=main)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ooalshei.github.io/RedCarD.jl/stable)
[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ooalshei.github.io/RedCarD.jl/dev)

**Red**uctive **Car**tan **D**ecompositions of Pauli operator algebras — and the
fixed-depth quantum circuits they produce.

Simulating $e^{-iHt}$ normally costs circuit depth that grows with $t$. When the Pauli
operators in $H$ close into a Lie algebra that admits a Cartan decomposition
$\mathfrak{g} = \mathfrak{k} \oplus \mathfrak{m}$, that trade disappears: writing

$$
H = K h K^\dagger, \quad K = \prod_j e^{i \theta_j k_j}, \quad ik_j \in \mathfrak{k}, \quad ih \in \mathfrak{h} \subseteq \mathfrak{m},
$$

with $\mathfrak{h}$ a maximal commuting subalgebra gives

$$
e^{-iHt} = K e^{-iht} K^\dagger
$$

where $e^{-iht}$ is a single layer of independent rotations. The depth is the same for
$t = 1$ as for $t = 10^6$.

RedCarD.jl builds that object: it constructs the algebra, splits it, finds the Cartan
subalgebra, and solves for the angles.

## Installation

Neither package is registered, so install the dependency first:

```julia
using Pkg
Pkg.add(url="https://github.com/ooalshei/SymplecticPauli.jl")
Pkg.add(url="https://github.com/ooalshei/RedCarD.jl")
```

## Example

```julia
using RedCarD

H = hamiltonian("TFIM", 4, [1.0, 0.5])                 # -J Σ XᵢXᵢ₊₁ - Jg Σ Zᵢ
decomposition = involutionlessdecomp(PauliList(collect(keys(H)), 4))

length(decomposition.g), length(decomposition.k)       # (28, 12)
tostring(decomposition.h)                              # ["-Z--", "--Z-", "Z---", "---Z"]

result = optimizer(H, decomposition.h, decomposition.k; convergence_tol=1e-8)
```

The optimizer reports its progress as it runs, starting from random angles:

```
Converged in 1230 iterations.
Final relative error: 8.685465962697732e-9
```

`result.H` is $h$, and `result.generators`/`result.angles` are the $k_j$ and $\theta_j$ that build
$K$. Conjugating back reproduces the original Hamiltonian to the convergence tolerance:

```julia
back = ad(result.H, decomposition.k, result.angles)
maximum(abs(get(back, p, 0.0im) - get(H, p, 0.0im)) for p in union(keys(H), keys(back)))
# 1.2e-8
```

## The reductive algorithm

`optimizer` varies all $\dim \mathfrak{k}$ angles against one scalar cost. Rotosolve is a local method, so
past about four qubits that stalls. The reductive algorithm of the
[paper](https://arxiv.org/abs/2512.06070) — the one this package is named for — takes the
problem apart instead: a Hamiltonian that commutes with a *multiplicative generating set* of
$\mathfrak{h}$ commutes with all of it, and the elements of $\mathfrak{k}$ split by which generator each one moves.

```julia
bstrings  = subalgred(decomposition.h)                      # h from a few generators
fragments = fragmentedsubspaces(decomposition.k, bstrings)  # k, split by which one it moves
cleangenerators!(fragments, bstrings)

result = reductive_optimizer(H, bstrings, fragments; convergence_tol=1e-8)
```

Each stage rotates $H$ into the centralizer of one more generator while varying only the
fragment of $\mathfrak{k}$ that acts on it. The stages are small and independent, and each one's cost is
the expectation of a *single* Pauli string — measurable on hardware, which is what makes the
quantum-assisted version of the algorithm possible.

TFIM at $J = 1$, $g = 0.5$, tolerance `1e-8`, best of three random starts, one-shot capped at
20000 sweeps:

| Qubits | $\dim \mathfrak{g}$ | $\dim \mathfrak{k}$ | one-shot sweeps | evals | time | error | reductive sweeps | evals | time | error |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 4 | 28 | 12 | 1780 | 64080 | 0.017 s | 9.1e-9 | 90 | 1320 | 0.0004 s | 1.0e-9 |
| 6 | 66 | 30 | 20000 † | 1800000 | 0.81 s | 1.6e-2 | 240 | 5580 | 0.0022 s | 9.3e-9 |
| 8 | 120 | 56 | 20000 † | 3360000 | 1.91 s | 1.2e-2 | 450 | 13920 | 0.0073 s | 6.7e-9 |
| 10 | 190 | 90 | 20000 † | 5400000 | 4.63 s | 4.1e-3 | 800 | 32460 | 0.026 s | 8.8e-9 |

† hit the cap without reaching the tolerance, from all three starts.

Both routes give the same circuit depth — the fragments partition the same $\mathfrak{k}$. Only one of
them gets there past a handful of qubits.

## What is in the box

| Step | Functions |
|:--|:--|
| Model Hamiltonians (Ising, XY, TFIM, TFXY, MFIM, HIM, Heisenberg, XXZ, Gross–Neveu) | `hamiltonian` |
| Dynamical Lie algebra closure | `dla` |
| Cartan decompositions, with an involution or without one | `cartandecomp`, `involutionlessdecomp` |
| Involutions | `evenoddx`, `evenoddy`, `evenoddz`, `typeIorII`, `typeIII` |
| Cartan subalgebra, its multiplicative generators, and the matching `k` fragments | `subalgfind`, `subalgred`, `fragmentedsubspaces`, `cleangenerators!` |
| Rotosolve optimizers, one-shot and reductive | `optimizer`, `reductive_optimizer` |

Pauli operators themselves — symplectic bit strings, `PauliList`, `PauliSentence`, `com`,
`ad`, `tostring`, `tomatrix` — come from
[SymplecticPauli.jl](https://github.com/ooalshei/SymplecticPauli.jl) and are re-exported, so
`using RedCarD` is enough.

## Threading

Algebra closures and the optimizer use every thread the session has; start Julia with
`julia -t auto`. A threaded closure returns the same algebra as a serial one, but discovers
its elements in a different order, so the generator ordering — and the particular set of
angles the optimizer lands on — can differ between runs. Every such solution is equally
valid.

## Documentation

Full documentation, a worked tutorial and the API reference are built with
[Documenter.jl](https://documenter.juliadocs.org):

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

then open `docs/build/index.html`. Every example in the documentation is a doctest and is
verified during that build.

## Reference

The algorithms implemented here follow

>*Fixed Depth Hamiltonian Simulation via Cartan Decomposition*,
>[Phys. Rev. Lett. 129, 070501](https://doi.org/10.1103/PhysRevLett.129.070501).  
>*RedCarD: A Quantum Assisted Algorithm for Fixed-Depth Unitary Synthesis via Cartan Decomposition*,
>[arXiv:2512.06070](https://arxiv.org/abs/2512.06070).
