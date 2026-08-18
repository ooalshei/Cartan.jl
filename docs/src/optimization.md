```@meta
CurrentModule = RedCarD
```

# Finding the circuit

A Cartan decomposition says a $K$ exists with $H = K h K^\dagger$. Finding it means solving
for the angles.

## The cost function

For an element $v$ of the Cartan subalgebra, we minimize

```math
f(\vec{\theta}) = \left\langle K(\vec{\theta}) v K(\vec{\theta})^\dagger, H \right\rangle .
```

Its minimum is attained exactly when $iK^\dagger H K$ lies in $\mathfrak{h}$, which is the
statement being solved for. The distinct coefficients matter — they are what stops different
Cartan elements from being traded off against each other — and `optimizer` builds them
internally from a sequence of mutually irrational values.

## Rotosolve

As a function of one angle with all others fixed, $f$ is a pure sinusoid,

```math
f(\theta_j) = c_0 + c_1\cos 2\theta_j + c_2 \sin 2\theta_j ,
```

so its minimizer is available in closed form. [`optimizer`](@ref) sweeps the generators,
setting each angle to its exact optimum in turn, and repeats until the weight of $H$ outside
$\mathfrak{h}$ falls below `convergence_tol`.

```jldoctest opt
julia> using RedCarD, Random

julia> H = hamiltonian("TFIM", 4, [1.0, 0.5]);

julia> decomposition = involutionlessdecomp(PauliList(collect(keys(H)), 4));

julia> Random.seed!(7);

julia> result = redirect_stdout(devnull) do
           optimizer(H, decomposition.h, decomposition.k; convergence_tol=1e-8)
       end;

julia> result.error < 1e-8, length(result.H)
(true, 4)
```

The optimizer reports its progress on stdout, which is silenced here; `redirect_stdout` is
not needed in normal use.

### Controlling the run

```julia
optimizer(
    H, decomposition.h, decomposition.k,
    initangles;              # defaults to π·rand(length(generators))
    maxiter = 5_000,         # 0 means no limit
    convergence_tol = 1e-8,  # relative Hilbert-Schmidt weight left outside 𝔥
    coeff_tol = 1e-12,       # drop coefficients smaller than this while rotating
    itertrack = true,        # add `iterations` and `calls` to the result
    timetrack = true,        # add `time`
)
```

`coeff_tol` is the one to reach for on large systems: intermediate sentences fill in with
terms of negligible weight, and pruning them keeps both memory and time in check at a cost
you control. Convergence is monitored every ten sweeps, so `maxiter` is best set to a
multiple of ten.

Rotosolve is a local method. If it stalls at an error far above `convergence_tol`, restart
from a different `initangles` rather than waiting.

## The reductive optimizer

[`reductive_optimizer`](@ref) solves the same problem in stages. Instead of asking for all
of $\mathfrak{h}$ at once it works through the generators returned by [`subalgred`](@ref):
first rotate $H$ so that it commutes with $b_1$, then with $b_2$ while keeping $b_1$, and so
on. Every stage only varies the fragment of $\mathfrak{k}$ that acts on that generator, so
each subproblem is much smaller than the full one.

```jldoctest opt
julia> bstrings = subalgred(decomposition.h);

julia> fragments = fragmentedsubspaces(decomposition.k, bstrings);

julia> cleangenerators!(fragments, bstrings);

julia> Random.seed!(11);

julia> reduced = redirect_stdout(devnull) do
           reductive_optimizer(H, bstrings, fragments; convergence_tol=1e-8)
       end;

julia> reduced.error < 1e-7
true

julia> length.(reduced.angles)          # one angle block per Cartan generator
3-element Vector{Int64}:
 6
 4
 2
```

`reduced.angles` is a vector of vectors, one per element of `bstrings`; concatenating them
in order gives the full circuit. This staging is what the quantum-assisted algorithm of the
[paper](https://arxiv.org/abs/2512.06070) needs, since each stage's cost is measurable on
hardware.

## Reading the result

Both optimizers return a `NamedTuple`:

| Field | Present | Meaning |
|:--|:--|:--|
| `H` | always | the transformed Hamiltonian, supported on $\mathfrak{h}$ |
| `generators` | always | the $\mathfrak{k}$ elements used, fragmented for the reductive form |
| `angles` | always | the solved angles, in the same shape as `generators` |
| `error` | always | relative Hilbert-Schmidt weight left outside $\mathfrak{h}$ |
| `iterations`, `calls` | `itertrack=true` | sweeps and cost evaluations |
| `time` | `timetrack=true` | seconds spent |

To rebuild $H$ from a result, conjugate back with the same generators and angles:

```jldoctest opt
julia> back = ad(result.H, decomposition.k, result.angles);

julia> maximum(
           abs(get(back, p, 0.0im) - get(H, p, 0.0im)) for
           p in union(keys(H), keys(back))
       ) < 1e-7
true
```

`ad` applies the generators from the last to the first, so undoing a rotation
sequence on its own means reversing the generators *and* the angles:
`ad(x, reverse(gens), -reverse(angles))`.
