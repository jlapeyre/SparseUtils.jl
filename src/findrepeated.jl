@inline _element_span(indhi, indlo) = indhi - indlo + 1

"""
    findrepeated(iter, min_seq_length::Integer)

Return an `Array` containing the indices of each maximal sequence of repeated elements in `iter` of length
not less than `min_seq_length`. The range of indices of each sequence is given by a `UnitRange`.
By construction, the output array of `UnitRanges` is sorted.

## Example
```julia
julia> a = [fill(1,3); fill("cat", 4); fill(0, 2);  5 ; fill(1,10)];

julia> findrepeated(a, 1)
5-element Array{Any,1}:
 1:3
 4:7
 8:9
 10:10
 11:20

julia> findrepeated(a, 2)
4-element Array{Any,1}:
 1:3
 4:7
 8:9
 11:20
```
"""
function findrepeated(iter, min_seq_length::Integer)
    inds = Vector{Any}()
    isempty(iter) && return inds
    newind = 1
    newx = first(iter)
    j = 0
    for (i, x) in enumerate(iter)
        j = i
        if x != newx
            last_seq_ind = i - 1
            if _element_span(last_seq_ind, newind) >= min_seq_length
                push!(inds, newind:last_seq_ind)
            end
            newind = i
            newx = x
        end
    end
    last_seq_ind = j
    if _element_span(last_seq_ind, newind) >= min_seq_length
        push!(inds, newind:last_seq_ind)
    end
    return inds
end

function compressinds!(inds)
    isempty(inds) && return inds
    inds[1] = inds[1] .- (inds[1].start - 1)
    for i in 2:length(inds)
        gap = inds[i].start - inds[i-1].stop - 1
        if gap > 0
            inds[i] = inds[i] .- gap
        end
    end
    return inds
end
compressinds(inds) = compressinds!(copy(inds))

function shiftrange(range, newstart::Integer=1)
    return range .- (range.start - newstart)
end


nothing
