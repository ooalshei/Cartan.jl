```@meta
CurrentModule = RedCarD
```

# Tutorial

This page runs the whole pipeline on one model — the transverse-field Ising chain on four
qubits — and then checks that the answer is what it claims to be.

```math
H = -J\sum_{i} X_i X_{i+1} - Jg \sum_i Z_i
```

## 1. The Hamiltonian

[`hamiltonian`](@ref) builds the standard models directly as a
`PauliSentence`: a dictionary from Pauli strings to coefficients.

```jldoctest tutorial
julia> using RedCarD

julia> H = hamiltonian("TFIM", 4, [1.0, 0.5]);

julia> length(H)
7

julia> sort(collect(keys(tostring(H))))
7-element Vector{String}:
 "---Z"
 "--XX"
 "--Z-"
 "-XX-"
 "-Z--"
 "XX--"
 "Z---"
```

Three `XX` bonds and four `Z` fields, as expected on an open chain. `-` marks an identity.

## 2. The algebra it generates

The Paulis in `H` do not commute, so repeatedly taking commutators generates more of them.
[`dla`](@ref) closes that process and returns the *dynamical Lie algebra*: the smallest set
of Pauli strings containing the generators and closed under commutation.

```jldoctest tutorial
julia> generators = PauliList(collect(keys(H)), 4);

julia> length(dla(generators))
28
```

Seven generators grow into a 28-dimensional algebra — and, crucially, they stop growing.
That is what makes a Cartan decomposition possible at all; an algebra that filled the whole
$4^Q$-dimensional space would give no advantage.

## 3. Splitting the algebra

A Cartan decomposition is a split $\mathfrak{g} = \mathfrak{k} \oplus \mathfrak{m}$ with

```math
[\mathfrak{k},\mathfrak{k}] \subseteq \mathfrak{k}, \qquad
[\mathfrak{k},\mathfrak{m}] \subseteq \mathfrak{m}, \qquad
[\mathfrak{m},\mathfrak{m}] \subseteq \mathfrak{k},
```

and it must put $iH$ in $\mathfrak{m}$. [`involutionlessdecomp`](@ref) finds one without
being told which involution to use: it grows the algebra from the generators while carrying
a sign, forced to $-1$ on everything it started from, and propagates that sign through every
commutator. If no string is ever reached with both signs, the assignment *is* a Cartan
decomposition.

```jldoctest tutorial
julia> decomposition = involutionlessdecomp(generators);

julia> length(decomposition.g), length(decomposition.k), length(decomposition.m)
(28, 12, 16)

julia> tostring(decomposition.h)
4-element Vector{String}:
 "-Z--"
 "--Z-"
 "Z---"
 "---Z"
```

`decomposition.h` is a Cartan subalgebra: a maximal commuting set inside $\mathfrak{m}$,
here the four single-qubit $Z$ operators. Rotating `H` into it is the whole game.

If the sign assignment had run into a contradiction, `k` and `m` would come back as
`nothing` — see [Algebras and decompositions](@ref) for what to do then.

## 4. Solving for the circuit, all at once

[`optimizer`](@ref) looks for the angles $\theta_j$ that rotate `H` into
`decomposition.h`, using [rotosolve](https://quantum-journal.org/papers/q-2021-01-28-391/):
each angle in turn is set to its exact optimum with the others held fixed, which is possible
because the cost is a pure sinusoid in every single angle.

```jldoctest tutorial
julia> using Random; Random.seed!(7);

julia> result = redirect_stdout(devnull) do          # the optimizer reports progress
           optimizer(H, decomposition.h, decomposition.k; convergence_tol=1e-8)
       end;

julia> result.error < 1e-8
true

julia> length(result.H)          # H is now supported on the 4 commuting elements
4
```

The returned `NamedTuple` carries everything the circuit needs:

| Field | Meaning |
|:--|:--|
| `result.H` | $h$, supported on the Cartan subalgebra |
| `result.generators` | the $k_j \in \mathfrak{k}$ |
| `result.angles` | the $\theta_j$ |
| `result.error` | relative Hilbert-Schmidt weight left outside $\mathfrak{h}$ |

That worked, but it took around 1800 sweeps over all twelve angles at once. Rotosolve is a
*local* method, and this cost function has all $\dim \mathfrak{k}$ angles coupled through a
single scalar. Four qubits is about where that stops being true: at six the same call runs
20000 sweeps and still sits three orders of magnitude short of the tolerance.

## 5. The reductive algorithm

The [paper](https://arxiv.org/abs/2512.06070) fixes this by refusing to ask for all of
$\mathfrak{h}$ at once. Two observations make that possible.

First, $\mathfrak{h}$ does not need $\dim \mathfrak{h}$ conditions imposed on it. Its elements
are products of a smaller *multiplicative* generating set $\{b_1,\dots,b_r\}$, and a
Hamiltonian that commutes with every $b_i$ commutes with all of $\mathfrak{h}$. For TFIM/TFXY,
$\mathfrak{h}$ will not be reduced. Beyond free-fermionic models, such as Heisenberg, we can drop
some terms. [`subalgred`](@ref) finds that set:

```jldoctest tutorial
julia> bstrings = subalgred(decomposition.h);

julia> tostring(bstrings)
4-element Vector{String}:
 "-Z--"
 "--Z-"
 "Z---"
 "---Z"
```

Second, the elements of $\mathfrak{k}$ split by *which* $b_i$ they first fail to commute
with — the pieces $\mathfrak{k}_i$ of the [home page](@ref RedCarD.jl).
[`fragmentedsubspaces`](@ref) performs that split, and
[`cleangenerators!`](@ref) drops any $b_i$ that no generator acts on:

```jldoctest tutorial
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

Twelve generators, but now in three groups: $\lvert\mathfrak{k}_1\rvert = 6$,
$\lvert\mathfrak{k}_2\rvert = 4$, $\lvert\mathfrak{k}_3\rvert = 2$.
[`reductive_optimizer`](@ref) walks through them: rotate $H$ until it commutes with $b_1$,
varying only the six generators that move $b_1$; then until it commutes with $b_2$, varying
only the four that move $b_2$ and leave $b_1$ alone; then $b_3$. Each stage is a small,
independent problem, and because $\mathfrak{k}_i$ commutes with $b_1,\ldots,b_{i-1}$, no
stage can undo the one before it.

```jldoctest tutorial
julia> Random.seed!(11);

julia> reduced = redirect_stdout(devnull) do
           reductive_optimizer(H, bstrings, fragments; convergence_tol=1e-8)
       end;

julia> reduced.error < 1e-7
true

julia> length.(reduced.angles)      # three small problems, not one twelve-angle problem
3-element Vector{Int64}:
 6
 4
 2

julia> length(reduced.H)            # same destination: H lands in 𝔥
4
```

There is a second reason to prefer this. Each stage minimizes
$f_r = \langle Kb_rK^\dagger, H_{r-1} \rangle$, the expectation value of a *single*
Pauli string, which is something a quantum computer can measure directly — that is what makes
the quantum-assisted version of the algorithm possible. The one-shot cost, a weighted sum
over all of $\mathfrak{h}$, is not.

## 6. How the two compare

Same model, same tolerance ($10^{-8}$), best of three random starts each, with the one-shot
optimizer capped at 20000 sweeps:

| Qubits | $\dim\mathfrak{g}$ | $\dim\mathfrak{k}$ | one-shot sweeps | cost evaluations | time | error | reductive sweeps | cost evaluations | time | error |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 4 | 28 | 12 | 1780 | 64080 | 0.017 s | 9.1e-9 | 90 | 1320 | 0.0004 s | 1.0e-9 |
| 6 | 66 | 30 | 20000<sup>†</sup> | 1800000 | 0.81 s | 1.6e-2 | 240 | 5580 | 0.0022 s | 9.3e-9 |
| 8 | 120 | 56 | 20000<sup>†</sup> | 3360000 | 1.91 s | 1.2e-2 | 450 | 13920 | 0.0073 s | 6.7e-9 |
| 10 | 190 | 90 | 20000<sup>†</sup> | 5400000 | 4.63 s | 4.1e-3 | 800 | 32460 | 0.026 s | 8.8e-9 |

<sup>†</sup> hit the iteration cap without reaching the tolerance, from all three starts.

The reductive algorithm converges at every size tested, and where both converge it needs
about 20× fewer sweeps and 50× fewer cost evaluations. The gap widens with system size,
because the one-shot problem couples more angles as $\dim \mathfrak{k}$ grows while the
reductive stages stay the same size — the largest fragment is $2(Q-1)$ generators regardless
of how big $\mathfrak{k}$ gets.

Reproduce it with:

```julia
using RedCarD, Random
for Q in (4, 6, 8, 10)
    H = hamiltonian("TFIM", Q, [1.0, 0.5], UInt32)
    d = involutionlessdecomp(PauliList(collect(keys(H)), Q))
    b = subalgred(d.h)
    frag = fragmentedsubspaces(d.k, b)
    cleangenerators!(frag, b)
    Random.seed!(1)
    full = optimizer(H, d.h, d.k; convergence_tol=1e-8, maxiter=20000,
                     itertrack=true, timetrack=true)
    Random.seed!(1)
    red = reductive_optimizer(H, b, frag; convergence_tol=1e-8, maxiter=20000,
                              itertrack=true, timetrack=true)
    @show Q full.iterations full.error red.iterations red.error
end
```

The cost-evaluation counts are the resource that matters on hardware: they assume the
three-point sampling rotosolve needs when the cost is measured rather than computed. In this
package the same information is obtained analytically in one pass, so the wall-clock ratio is
smaller than the evaluation ratio.

## 7. Checking it

Nothing so far proves the factorization is correct, so verify it directly. Conjugating the
result back with the same generators and angles has to return the original Hamiltonian:

```jldoctest tutorial
julia> back = ad(result.H, decomposition.k, result.angles);   # K h K†

julia> maximum(
           abs(get(back, p, 0.0im) - get(H, p, 0.0im)) for
           p in union(keys(H), keys(back))
       ) < 1e-7
true
```

The reductive result factorizes the same way once its per-stage generators and angles are
concatenated in order:

```jldoctest tutorial
julia> generators = vcat(reduced.generators...); angles = vcat(reduced.angles...);

julia> backr = ad(reduced.H, generators, angles);

julia> maximum(
           abs(get(backr, p, 0.0im) - get(H, p, 0.0im)) for
           p in union(keys(H), keys(backr))
       ) < 1e-7
true
```

And since $H$ and $h$ are related by a unitary, they must have the same spectrum:

```jldoctest tutorial
julia> using LinearAlgebra

julia> isapprox(
           sort(eigvals(Hermitian(tomatrix(H)))),
           sort(eigvals(Hermitian(tomatrix(result.H))));
           atol=1e-6,
       )
true

julia> round(minimum(eigvals(Hermitian(tomatrix(H)))), digits=8)
-3.42703409
```

The ground-state energy is now available from `result.H` alone: it is a sum of four
commuting terms, so its extremal eigenvalues are just $\pm$ the sum of the absolute values
of its coefficients.

## 8. Why this is a fixed-depth circuit

`result` describes

```math
e^{-iHt} = K e^{-iht} K^\dagger,
\qquad K = \prod_{j=1}^{12} e^{i\theta_j k_j},
```

so the circuit is 12 rotations, 4 commuting rotations whose angles are the only thing that
depends on $t$, and 12 rotations back — from either optimizer, since the reductive stages
partition the same $\mathfrak{k}$:

```jldoctest tutorial
julia> length(result.angles), length(result.H)
(12, 4)

julia> sum(length, reduced.angles), length(reduced.H)
(12, 4)
```

Evolving to $t = 10^{6}$ costs exactly as much as evolving to $t = 1$. What the reductive
algorithm buys is not a shorter circuit — it is being able to find one at all past a handful
of qubits.

## Where to go next

* [Algebras and decompositions](@ref) — involutions, the involutionless search, and what to
  do when no decomposition exists.
* [Finding the circuit](@ref) — the cost function both optimizers minimize, and the knobs
  that matter on larger systems.
