"""
    module IJV

Algorithms for sparse matrices stored in COO form.

`IJV` does not depend on any types not in `Core`.
Do not import any symbol from `IJV` into another module. The methods
in `IJV` are incompatible with functions in `Base` with the same
(unqualified) name. Instead, use `IJV` like this

#Example

```julia
Base.getindex(S::MySparseMatrixCOO, i::Integer, j::Integer) = IJV.getindex(S.I, S.J, S.V, i, j)
```
"""
module IJV

import Random
import StatsBase
import SparseArrays
import SparseUtils.findrepeated

spzeros(m::Integer, n::Integer) = spzeros(Float64, m, n)
spzeros(::Type{Tv}, m::Integer, n::Integer) where {Tv} = spzeros(Tv, Int, m, n)
function spzeros(::Type{Tv}, ::Type{Ti}, m::Integer, n::Integer) where {Tv, Ti}
    ((m < 0) || (n < 0)) && throw(ArgumentError("invalid Array dimensions"))
    return (m, n, Ti[], Ti[], Tv[])
#    SparseMatrixCSC(m, n, fill(one(Ti), n+1), Vector{Ti}(), Vector{Tv}())
end

### sorting and indexing

"""
    issorted(I::AbstractVector, J::AbstractVector)

Return true if the COO indices `I` and `J` are sorted canonically.
"""
issorted(I::AbstractVector, J::AbstractVector) = Base.issorted(zip(J, I))

function getindex(I, J, V, i::Integer, j::Integer; zeroel=zero(eltype(V)))
    # colrange = searchsorted(J, j)
    # isempty(colrange) && return zeroel
    # rowrange = searchsorted(view(I, colrange), i)
    # isempty(rowrange) && return zeroel
    # index = first(rowrange) + first(colrange) - 1
    index = findindex(I, J, i, j)
    index == nothing && return zeroel
    return V[index]
end

@inline function findindex(I::AbstractArray, J::AbstractArray, i::Integer, j::Integer)
    colrange = searchsorted(J, j)
    isempty(colrange) && return nothing
    rowrange = searchsorted(view(I, colrange), i)
    isempty(rowrange) && return nothing
    return first(rowrange) + first(colrange) - 1
end

function ijvinsert!(ind, arrays = (I, J, V), elements = (i, j, val))
    for (array, element) in zip(arrays, elements)
        insert!(array, ind, element)
    end
    return nothing
end

function ijvsetindex!(ind, arrays = (I, J, V), elements = (i, j, val))
    for (array, element) in zip(arrays, elements)
        array[ind] = element
    end
    return nothing
end

function setindex!(I, J, V, val, i::Integer, j::Integer; zeroel=zero(eltype(V)))
    rowrange = searchsorted(I, i)
    if isempty(rowrange)
        ind = first(rowrange)
        ijvinsert!(ind, (I, J, V), (i, j, val))
    else
        colrange = searchsorted(view(J, rowrange), j)
        if isempty(colrange)
            ijvinsert!(first(colrange), (I, J, V), (i, j, val))
        else
            valind = first(colrange) + first(rowrange) - 1
            ijvsetindex!(valind, (I, J, V), (i, j, val))
        end
    end
    return val
end

sortpermindex(I, J) = sortperm(1:length(I), by = i -> @inbounds (J[i], I[i]))
function ijvsort!(I, J, V)
    p = sortpermindex(I, J)
    permute!(I, p)
    permute!(J, p)
    permute!(V, p)
    return nothing
end

function permutedims!(m, n, I, J, V)
    ijvsort!(J, I, V) # swap I and J before sorting
    return (n, m, J, I, V)
end
permutedims(m, n, I, J, V) = permutedims!(m, n, copy(I), copy(J), copy(V))

function Array(m, n, I, J, V)
    a = zeros(eltype(V), m, n)
    for ind in 1:length(I)
        a[I[ind], J[ind]] = V[ind]
    end
    return a
end

### dropzeros, dropstored!

function dropzeros!(I, J, V::AbstractArray{T}) where T
    mask = V .== zero(T)
    map(x -> deleteat!(x, mask), (I, J, V))
end
dropzeros(I, J, V::AbstractArray) = dropzeros!(copy(I), copy(J), copy(V))

function dropstored!(I, J, V, i::Integer, j::Integer)
    index = findindex(I, J, i, j)
    index == nothing && return nothing
    deleteat!(I, index)
    deleteat!(J, index)
    deleteat!(V, index)
    return nothing
end

### prunecols!

# function prunecols!(J, min_entries)
#     col_entry_counts = StatsBase.countmap(J)
#     keepcols = filter(k -> col_entry_counts[k] >= min_entries, keys(col_entry_counts))
#     old_col_nums = collect(keepcols)
#     sort!(old_col_nums) # should not be necessary. We can avoid allocating if key avoid sort
#     renumber_dict = Dict(col_num => i for (i, col_num) in enumerate(old_col_nums))
#     @show renumber_dict
#     for i in 1:length(J)
# #        @show i
#         J[i] = renumber_dict[J[i]]
#     end
#     return J
# end

## We'd like `x[inds]`, but this appears to be a missing feature in Iterators.flatten
_applyinds(x, inds) = [x[i] for i in inds]

function prunecols(m, n, I, J, V, min_entries)
    inds = Iterators.flatten(findrepeated(J, min_entries))
    return (m, n, map(x -> _applyinds(x, inds), (I, J, V))...)
end

### rotation

# This is copied from sparsematrix.jl, modified only to change the order of the arguments
function rot180!(m, n, I, J, V)
   @inbounds for i=1:length(I)
        I[i] = m - I[i] + 1
        J[i] = n - J[i] + 1
    end
    return (m, n, I, J, V)
end

function rotr90!(m, n, I, J, V)
    #old col inds are new row inds
   @inbounds for i=1:length(I)
        I[i] = m - I[i] + 1
    end
    return (n, m, J, I, V)
end

function rotl90!(m, n, I, J, V)
    #old row inds are new col inds
  @inbounds for i=1:length(J)
        J[i] = n - J[i] + 1
    end
    return (n, m, J, I, V)
end

### sprand

_empty_IJV(T) = (Int[], Int[], T[])

function sprand(m::Integer, n::Integer, density::AbstractFloat,
                rfn::Function, ::Type{T}=eltype(rfn(rng,1))) where T
    return sprand(Random.GLOBAL_RNG, m, n, density, rfn, eltype)
end

function sprand(rng::Random.AbstractRNG, m::Integer, n::Integer, density::AbstractFloat,
                rfn::Function, ::Type{T}=eltype(rfn(rng,1))) where T
    N = m*n
    if N == 0
        (I, J, V) = _empty_IJV(T)
    elseif N == 1
        if rand(rng) <= density
            (I, J, V) = ([1], [1], rfn(rng, 1))
        else
            (I, J, V) = _empty_IJV(T)
        end
    else
        I, J = SparseArrays.sprand_IJ(rng, m, n, density)
        V = rfn(rng, length(I))
    end
    return (m, n, I, J, V)
end

truebools(r::Random.AbstractRNG, n::Integer) = fill(true, n)
sprand(m::Integer, n::Integer, density::AbstractFloat) = sprand(Random.GLOBAL_RNG,m,n,density)

sprand(r::Random.AbstractRNG, m::Integer, n::Integer, density::AbstractFloat) = sprand(r,m,n,density,rand,Float64)
sprand(r::Random.AbstractRNG, ::Type{T}, m::Integer, n::Integer, density::AbstractFloat) where {T} = sprand(r,m,n,density,(r, i) -> rand(r, T, i), T)
sprand(r::Random.AbstractRNG, ::Type{Bool}, m::Integer, n::Integer, density::AbstractFloat) = sprand(r,m,n,density, truebools, Bool)
sprand(::Type{T}, m::Integer, n::Integer, density::AbstractFloat) where {T} = sprand(Random.GLOBAL_RNG, T, m, n, density)

end # module
