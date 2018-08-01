"""
    countmap(iter; datatype=Any, dicttype=Dict)

Return an `AbstractDict` counting the number of occurences of elements
in `iter`.

Unlike `StatsBase.countmap`, `DataUtils.countmap`
does not require that one first calls `collect` on `iter` if
`iter` is not an `AbstractArray`. `DataUtils.countmap` accepts
any iterator, including one that may generate more values than can
be stored in memory.

`DataUtils.countmap` does not implement a fast, radix sort
as `StatsBase.countmap` does.

# Example
```julia
julia> import StatsBase; import DataUtils;

julia> d = Dict(rand(1:1000) => rand(1:5) for i in 1:100);

julia> StatsBase.countmap(values(d))
ERROR: MethodError: no method matching countmap(::Base.ValueIterator{Dict{Int64,Int64}})
...

julia> DataUtils.countmap(values(d))
Dict{Any,Int64} with 5 entries:
  4 => 19
  2 => 21
  3 => 20
  5 => 11
  1 => 25
```
"""
function countmap(iter; datatype=Any, dicttype=Dict)
    return countmap!(dicttype{datatype,Int}(), iter)
end

"""
    countmap(iter::AbstractArray{T}; dicttype=Dict) where {T}

Equivalent to calling `countmap(iter; datatype=Any, dicttype=Dict)`,
with `datatype=T`.
"""
function countmap(a::AbstractArray{T}; dicttype=Dict) where {T}
    return countmap!(dicttype{T,Int}(), a)
end


function countmap!(cm::Dict, iter)
    for v in iter
        index = Base.ht_keyindex2!(cm, v)
        if index > 0
            @inbounds cm.vals[index] += 1
        else
            @inbounds Base._setindex!(cm, 1, v, -index)
        end
    end
    return cm
end

"""
    countmap!(cm::AbstractDict, iter)

Add counts of the values in `iter`
to the countmap `cm`. See `countmap`
"""
function countmap!(cm::AbstractDict, iter)
    for v in iter
        if haskey(cm, v)
            @inbounds cm[v] += 1
        else
            @inbounds cm[v] = 1
        end
    end
    return cm
end
