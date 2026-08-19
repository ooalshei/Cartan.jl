```@meta
CurrentModule = RedCarD
```

# RedCarD.jl

**Red**uctive **Car**tan **D**ecompositions of Pauli operator algebras — and the
fixed-depth quantum circuits they produce.

Simulating $e^{-iHt}$ normally costs circuit depth that grows with $t$. If the Pauli
operators in $H$ close into a Lie algebra $\mathfrak{g}$ that admits a Cartan decomposition
$\mathfrak{g} = \mathfrak{k} \oplus \mathfrak{m}$, that trade disappears. Writing

```math
H = K\,h\,K^\dagger, \qquad
K = \prod_j e^{i\theta_j \tilde{k}_j},\quad i\tilde{k}_j \in \mathfrak{k}, \quad ih \in \mathfrak{h},
```

with $\mathfrak{h} \subseteq \mathfrak{m}$ a *Cartan subalgebra* — a maximal set of mutually
commuting elements — the time evolution becomes

```math
e^{-iHt} = K\,e^{-iht}\,K^\dagger ,
```

and since everything inside $\mathfrak{h}$ commutes, $e^{-iht}$ is a layer of independent
rotations. The depth is set by $\dim \mathfrak{k}$ and $\dim \mathfrak{h}$,
and is the same for $t = 1$ as for $t = 10^6$.

RedCarD.jl builds that object: it constructs the algebra, splits it, finds the Cartan
subalgebra, and solves for the angles $\theta_j$.

## Why *reductive*

The existence of $K$ is a theorem; finding it is an optimization. The standard route fixes an
element $v \in \mathfrak{h}$ whose exponential map $e^{sv}$ for $s \in \mathbb{R}$ is *dense* in
$\exp{\mathfrak{h}}$ and minimizes

```math
f(\vec{\theta})
  = i\left\langle K(\vec{\theta})vK(\vec{\theta})^\dagger, H \right\rangle,
\qquad
\langle A, B\rangle = \frac{1}{2^{Q}}\,\mathrm{Tr}\!\left[A B\right],
```

over all $\dim \mathfrak{k}$ angles at once. Its minimizer is the $K$ we want, but
every angle is coupled to every other through a single scalar, and the local optimizer used
to solve it stalls once $\mathfrak{k}$ grows past a few dozen elements.

The **reductive** decomposition — the one this package is named for — exploits something a
generic Lie algebra does not have: Pauli strings *multiply*. Every element of $\mathfrak{h}$
is a product of a small generating set

```math
\mathfrak{h} \subseteq \big\langle\, b_1, b_2, \ldots, b_r \,\big\rangle ,
\qquad ib_j \in \mathfrak{h} ,
```

with $r$ as small as $\log_2\!\big(\dim \mathfrak{h}+ 1\big)$, since $r$
independent strings already generate $2^r - 1$ of them. A Hamiltonian that commutes with
every $b_j$ therefore commutes with all of $\mathfrak{h}$, which replaces
$\dim \mathfrak{h}$ conditions by $r$ of them. The same generators split
$\mathfrak{k}$: each $ik \in \mathfrak{k}$ belongs to the piece indexed by the *first*
generator it fails to commute with,

```math
\mathfrak{k}_i = \big\{\, k \in \mathfrak{k} \;:\;
  [k, b_i] \neq 0 \;\text{ and }\; [k, b_j] = 0 \ \ \forall\, j < i \,\big\} ,
\qquad
\mathfrak{k} \supseteq \mathfrak{k}_1 \sqcup \cdots \sqcup \mathfrak{k}_r ,
```

the leftovers being the elements that commute with all of $\mathfrak{h}$. The circuit
factorizes to match,

```math
K = K_1 K_2 \cdots K_r , \qquad
K_j = \prod_{ik \in \mathfrak{k}_j} e^{i\alpha_k k} ,
```

and the angles are solved for one stage at a time. Stage $i$ minimizes

```math
f_r(\vec{\alpha})
  = \left\langle K_r(\vec{\alpha}) b_r K_r(\vec{\alpha})^\dagger, H_{r-1}\right\rangle ,
\qquad H_0 = H ,
```

and hands on $H_r$, the part of the rotated Hamiltonian that commutes with
$b_1,\ldots,b_r$. Because $\mathfrak{k}_r$ commutes with every earlier generator, no stage
can undo the one before it — the constraints accumulate.

Two things follow. The optimization is now $r$ small independent problems instead of one
large coupled one, which is what makes it converge at sizes where the dense cost function
does not. And each $f_r$ is the expectation value of a *single* Pauli string rather than a
weighted sum over all of $\mathfrak{h}$ — something a quantum computer can measure directly,
which is what the quantum-assisted form of the algorithm rests on.

## Installation

Neither package is registered, so install the dependency first:

```julia
using Pkg
Pkg.add(url="https://github.com/ooalshei/SymplecticPauli.jl")
Pkg.add(url="https://github.com/ooalshei/RedCarD.jl")
```

## Quick start

```jldoctest quickstart
julia> using RedCarD

julia> H = hamiltonian("TFIM", 4, [1.0, 0.5]);   # -J Σ XᵢXᵢ₊₁ - Jg Σ Zᵢ

julia> algebra = PauliList(collect(keys(H)), 4);

julia> decomposition = involutionlessdecomp(algebra);

julia> length(decomposition.g), length(decomposition.k), length(decomposition.m)
(28, 12, 16)

julia> tostring(decomposition.h)      # the commuting subalgebra H is rotated into
4-element Vector{String}:
 "-Z--"
 "--Z-"
 "Z---"
 "---Z"
```

## Solving for the angles

One call turns that decomposition into a circuit by varying every generator of
$\mathfrak{k}$ at once:

```julia
result = optimizer(H, decomposition.h, decomposition.k; convergence_tol = 1e-8)
```

```
Converged in 1230 iterations.
Final relative error: 8.685465962697732e-9
```

`result.H` is $h$, `result.generators` and `result.angles` are the $k_j$ and $\theta_j$
that build $K$.

The reductive route builds the $b_i$ and the $\mathfrak{k}_i$ of the previous section and
solves stage by stage:

```julia
bstrings  = subalgred(decomposition.h)                      # the generators bᵢ of 𝔥
fragments = fragmentedsubspaces(decomposition.k, bstrings)  # the pieces 𝔨ᵢ of 𝔨
cleangenerators!(fragments, bstrings)                       # drop any bᵢ nothing acts on

result = reductive_optimizer(H, bstrings, fragments; convergence_tol = 1e-8)
```

`result.generators` and `result.angles` come back fragmented the same way, one block per
$b_i$; concatenating them in order gives the same $K$.

Same model, same tolerance, best of three random starts (one-shot capped at 20000 sweeps):

| Qubits | $\dim \mathfrak{k}$ | one-shot sweeps | error | reductive sweeps | error |
|--:|--:|--:|--:|--:|--:|
| 4 | 12 | 1780 | 9.1e-9 | 90 | 1.0e-9 |
| 6 | 30 | 20000$^\dagger$ | 1.6e-2 | 240 | 9.3e-9 |
| 8 | 56 | 20000$^\dagger$ | 1.2e-2 | 450 | 6.7e-9 |
| 10 | 90 | 20000$^\dagger$ | 4.1e-3 | 800 | 8.8e-9 |

$^\dagger$ hit the cap without reaching the tolerance. Both routes produce the same circuit depth;
only one of them gets there. The [Tutorial](@ref) walks through both and checks that the
factorization really does reproduce `H`.

## The pipeline

| Step | What it does | Functions |
|:--|:--|:--|
| Build | model Hamiltonians as Pauli sentences | [`hamiltonian`](@ref) |
| Close | the Lie algebra generated by a set of Paulis | [`dla`](@ref) |
| Split | $\mathfrak{g} = \mathfrak{k} \oplus \mathfrak{m}$ | [`cartandecomp`](@ref), [`involutionlessdecomp`](@ref) |
| Reduce | pick $\mathfrak{h}$, its multiplicative generators, and the matching $\mathfrak{k}$ fragments | [`subalgfind`](@ref), [`subalgred`](@ref), [`fragmentedsubspaces`](@ref), [`cleangenerators!`](@ref) |
| Solve | for the circuit angles, all at once or one Cartan generator at a time | [`optimizer`](@ref), [`reductive_optimizer`](@ref) |

## Where things live

Pauli operators are symplectic bit strings: `2Q` bits per operator, the low half marking
`X` support and the high half `Z` support. That representation, the containers built on it
(`PauliList`, `PauliSentence`) and the operations on them (`com`,
`ad`, `tostring`, `tomatrix`) come from
[SymplecticPauli.jl](https://github.com/ooalshei/SymplecticPauli.jl) and are re-exported
here, so `using RedCarD` is all you need.

## Threading

The algebra closures and the optimizer parallelize themselves when the session has more than
one thread; start Julia with `julia -t auto` to use it. The algebra a closure returns is the
same either way, but its elements are discovered in a different order, so the generator
ordering — and with it the particular set of angles the optimizer lands on — can differ
between a serial and a threaded run. Every such solution is equally valid.

## Reference

The algorithms implemented here follow

>*Fixed Depth Hamiltonian Simulation via Cartan Decomposition*,
>[Phys. Rev. Lett. 129, 070501](https://doi.org/10.1103/PhysRevLett.129.070501).  
>*RedCarD: A Quantum Assisted Algorithm for Fixed-Depth Unitary Synthesis via Cartan Decomposition*,
>[arXiv:2512.06070](https://arxiv.org/abs/2512.06070).
