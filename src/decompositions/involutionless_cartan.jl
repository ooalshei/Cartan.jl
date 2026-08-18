function _dla(string::Pair{<:Unsigned,<:Integer}, iter, Q::Integer)
    strings = Pair.(com.(string.first, first.(iter), Q), string.second .* last.(iter))
    return filter!(x -> !x.first.second, strings) .|> (x -> Pair(x.first.first, x.second))
end

"""
    _BitSigns

Dense `string => ±1` map backed by two bitmaps, the counterpart of [`_newmembers`](@ref) for
the signed closure: `seen` records whether a string has been reached at all and `positive`
which of the two subspaces it was reached in.
"""
struct _BitSigns
    seen::BitVector
    positive::BitVector
end
_BitSigns(n::Integer) = _BitSigns(falses(n), falses(n))

_newsigns(::Type{T}, Q::Integer) where {T<:Unsigned} =
    2 * Q <= _maxbitmap ? _BitSigns(1 << (2 * Q)) : Dict{T,Int8}()

@inline function _signof(signs::_BitSigns, string::Unsigned)
    i = Int(string) + 1
    @inbounds signs.seen[i] || return Int8(0)
    @inbounds return signs.positive[i] ? Int8(1) : Int8(-1)
end
@inline function _setsign!(signs::_BitSigns, string::Unsigned, sign::Int8)
    i = Int(string) + 1
    @inbounds signs.seen[i] = true
    @inbounds signs.positive[i] = sign > 0
    return sign
end
@inline _signof(signs::AbstractDict, string::Unsigned) = get(signs, string, Int8(0))
@inline _setsign!(signs::AbstractDict, string::Unsigned, sign::Int8) =
    (signs[string] = sign)

"""
    _newsignedalgebra(paulis, k)

Start a signed closure from `paulis`: every generator is tagged `-1` (it belongs to
``\\mathfrak{m}``) unless it appears in `k`, in which case it is tagged `+1`.

Returns the tagged algebra together with the `string => sign` map used to detect both
repeats and sign contradictions. Keeping that map is what lets the closure avoid the
`unique!` rescan of the whole algebra after every commutator batch.
"""
function _newsignedalgebra(paulis::PauliList{T,Q}, k) where {T,Q}
    algebra = Pair{T,Int8}[]
    sizehint!(algebra, length(paulis))
    signs = _newsigns(T, Q)
    tag! =
        (p, sign) ->
            iszero(_signof(signs, p)) &&
            (_setsign!(signs, p, sign); push!(algebra, p => sign))
    if isnothing(k)
        for p in paulis
            tag!(p, Int8(-1))
        end
    else
        inner = Set{T}(k)
        for p in paulis
            p in inner && tag!(p, Int8(1))
        end
        for p in paulis
            p in inner || tag!(p, Int8(-1))
        end
    end
    return algebra, signs
end

"""
    _absorbsigned!(algebra, signs, newpaulis) -> contradiction

Append the `string => sign` pairs of `newpaulis` that are new to `algebra`, and report
whether any of them repeats a string already carrying the opposite sign. Such a repeat means
no consistent involution exists, so the first sign seen is kept and the caller falls back to
closing the algebra without signs.
"""
@inline function _tagsigned!(
    algebra::Vector{Pair{T,Int8}},
    signs,
    string::T,
    sign::Int8,
) where {T}
    previous = _signof(signs, string)
    if iszero(previous)
        _setsign!(signs, string, sign)
        push!(algebra, string => sign)
        return false
    end
    return previous != sign  # same string in both subspaces: no involution can exist
end

function _absorbsigned!(algebra::Vector{Pair{T,Int8}}, signs, newpaulis) where {T}
    contradiction = false
    for pair in newpaulis
        contradiction |= _tagsigned!(algebra, signs, T(pair.first), Int8(pair.second))
    end
    return contradiction
end

function _signedcommutators!(
    algebra::Vector{Pair{T,Int8}},
    signs,
    pair::Pair{T,Int8},
    indices::AbstractRange,
    Q::Integer,
) where {T}
    contradiction = false
    for j in indices
        other = @inbounds algebra[j]
        product = com(pair.first, other.first, Q)
        product.second && continue
        contradiction |=
            _tagsigned!(algebra, signs, T(product.first), Int8(pair.second * other.second))
    end
    return contradiction
end

function _signedresult(
    algebra::Vector{Pair{T,Int8}},
    Q::Integer,
    contradiction::Bool,
) where {T}
    alg = PauliList{T,Q}(first.(algebra), iscopy=false, check=false)
    contradiction && return (g=alg, h=subalgfind(alg), k=nothing, m=nothing)
    m = PauliList{T,Q}(
        [p.first for p in algebra if p.second == -1],
        iscopy=false,
        check=false,
    )
    k = PauliList{T,Q}(
        [p.first for p in algebra if p.second == 1],
        iscopy=false,
        check=false,
    )
    return (g=alg, k=k, m=m, h=subalgfind(m))
end

function _serial_involutionlessdecomp(
    paulis::PauliList{T,Q},
    startind::Integer=1;
    k=nothing,
) where {T,Q}
    signedalgebra, signs = _newsignedalgebra(paulis, k)
    contradiction, converged, initialind =
        _serial_signedclosure!(signedalgebra, signs, Int(startind), Q)
    converged && return _signedresult(signedalgebra, Q, contradiction)

    # A contradiction rules out the involutionless splitting, so only the algebra itself is
    # closed from here on.
    algebra, members =
        _newalgebra(PauliList{T,Q}(first.(signedalgebra), iscopy=false, check=false))
    _serial_closure!(algebra, members, initialind)
    return (g=algebra, h=subalgfind(algebra), k=nothing, m=nothing)
end

function _serial_signedclosure!(
    signedalgebra::Vector{Pair{T,Int8}},
    signs,
    initialind::Int,
    Q::Integer,
) where {T}
    finalind = length(signedalgebra)
    contradiction = false
    while !contradiction
        for i in finalind:-1:initialind
            contradiction |=
                _signedcommutators!(signedalgebra, signs, signedalgebra[i], (i-1):-1:1, Q)
        end
        finalind == length(signedalgebra) && return contradiction, true, initialind
        initialind = finalind + 1
        finalind = length(signedalgebra)
    end
    return contradiction, false, initialind
end

function _parallel_involutionlessdecomp(
    paulis::PauliList{T,Q},
    startind::Integer=1;
    k=nothing,
) where {T,Q}
    signedalgebra, signs = _newsignedalgebra(paulis, k)
    contradiction, converged, initialind =
        _parallel_signedclosure!(signedalgebra, signs, Int(startind), Q)
    converged && return _signedresult(signedalgebra, Q, contradiction)

    algebra, members =
        _newalgebra(PauliList{T,Q}(first.(signedalgebra), iscopy=false, check=false))
    _parallel_closure!(algebra, members, initialind)
    return (g=algebra, h=subalgfind(algebra), k=nothing, m=nothing)
end

function _parallel_signedclosure!(
    signedalgebra::Vector{Pair{T,Int8}},
    signs,
    initialind::Int,
    Q::Integer,
) where {T}
    finalind = length(signedalgebra)
    contradiction = false
    while !contradiction
        for i in finalind:-1:initialind
            if i - 1 < _minchunk
                contradiction |= _signedcommutators!(
                    signedalgebra,
                    signs,
                    signedalgebra[i],
                    (i-1):-1:1,
                    Q,
                )
            else
                iter = view(signedalgebra, (i-1):-1:1)
                chunks =
                    Iterators.partition(iter, max(1, length(iter) ÷ Threads.nthreads()))
                tasks = map(chunks) do chunk
                    Threads.@spawn _dla(signedalgebra[i], chunk, Q)
                end
                for task in tasks
                    contradiction |= _absorbsigned!(signedalgebra, signs, fetch(task))
                end
            end
        end
        finalind == length(signedalgebra) && return contradiction, true, initialind
        initialind = finalind + 1
        finalind = length(signedalgebra)
    end
    return contradiction, false, initialind
end

@doc raw"""
    involutionlessdecomp(paulis::PauliList)

Compute an involutionless Cartan decomposition of the algebra generated by `paulis`.

Build the dynamical Lie algebra ``\mathfrak{g}`` generated by `paulis` while tracking a sign
bit for involution consistency. Return a `NamedTuple` containing the generated algebra and
the subspaces ``\mathfrak{k}`` and ``\mathfrak{m}`` associated with the Cartan decomposition
that keeps `paulis` in ``\mathfrak{m}``. Return a Cartan subalgebra ``\mathfrak{h} \subseteq
\mathfrak{m}``. If no such decomposition exists, `k` and `m` are `nothing` and the returned
``\mathfrak{h}`` is a subalgebra of ``\mathfrak{g}``.

See also [`dla`](@ref), [`cartandecomp`](@ref), [`subalgfind`](@ref).

# Examples
```jldoctest
julia> ham = hamiltonian("TFIM", 4, [1, 1], UInt32);

julia> decomp = involutionlessdecomp(PauliList(collect(keys(ham)), 4));

julia> length(decomp.g), length(decomp.k), length(decomp.m), length(decomp.h)
(28, 12, 16, 4)

julia> println(decomp.h)
["-Z--", "--Z-", "Z---", "---Z"]
```
"""
involutionlessdecomp(paulis::PauliList) =
    if isone(Threads.nthreads())
        _serial_involutionlessdecomp(paulis, 1)
    else
        _parallel_involutionlessdecomp(paulis, 1)
    end
