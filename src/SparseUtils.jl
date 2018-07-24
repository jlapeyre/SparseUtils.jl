module SparseUtils

export c_to_julia_index, c_to_julia_index!, sparsity,
       transpose_concrete, numconnections, sparse_stats

import SparseArrays
import SparseArrays: SparseMatrixCSC, nnz
import DataStructuresUtils
import StatsBase
import Printf

"""
    c_to_julia_index!(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC

Convert a sparse matrix with zero-based indices to a `SparseMatrixCSC`.
`colptr0` and `rowval0` are altered in place. `nzval` is not copied.
"""
function c_to_julia_index!(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC
    colptr = colptr0 .+ 1
    rowval = rowval0 .+ 1
    m = length(colptr)-1
    n = maximum(rowval)
    return SparseArrays.SparseMatrixCSC(n, m, colptr, rowval, nzval)
end

"""
    c_to_julia_index(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC

Convert a sparse matrix with zero-based indices to a `SparseMatrixCSC`.
`colptr0` and `rowval0` are copied and altered. `nzval` is not copied.
"""
function c_to_julia_index(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC
    c_to_julia_index!(deepcopy(colptr0), deepcopy(rowval0), nzval)
end

"""
    sparsity(sp::SparseMatrixCSC)

Compute the sparisty of `sp`. This assumes that there
are no stored values equal to zero.
"""
function sparsity(sp::SparseMatrixCSC)
    return SparseArrays.nnz(sp) / length(sp)
end

"""
    transpose_concrete(M::SparseMatrixCSC)

Return the concrete transpose of the sparse matrix `M`. That is,
return a new `SparseMatrixCSC`, rather than `M` wrapped in `Transpose`.

Note that this does *not* return (or have any relation to) a dense matrix.
"""
function transpose_concrete(M::SparseMatrixCSC)
    I1, J, V = SparseArrays.findnz(M)
    m, n = size(M)
    return SparseArrays.sparse(J, I1, V, n, m)
end

"""
    nnz(sp::SparseMatrixCSC, colnum::Integer)

Return the number of structurally non-zero elements in the `colnum`th
column of `sp`.
"""
function SparseArrays.nnz(sp::SparseMatrixCSC, colnum::Integer)
    return length(SparseArrays.nzrange(sp, colnum))
end

"""
    numconnections(sp::SparseMatrixCSC; rev=true, byvalue=true)

Return a sorted countmap of the number of structural non-zeros
in the columns of `sp`. The keys controlling the sorting are
passed to `StatsBase.sort`.

# Examples
In this example, there are 132781 columns with one non-zero element,
32096 columns with two, etc.
```julia
julia> conns = numconnections(sp);

julia> collect(Iterators.take(conns,3))
3-element Array{Pair{Int64,Int64},1}:
 1 => 132781
 2 => 32096
 3 => 17478
```
"""
function numconnections(sp::SparseMatrixCSC; rev=true, byvalue=true)
    iter = (nnz(sp, i) for i in 1:size(sp)[2])
    spcountmap = DataStructuresUtils.countmap(iter, datatype=Int)
    sorted_spcountmap = StatsBase.sort(spcountmap, rev=rev, byvalue=rev)
    return sorted_spcountmap
end

"""
    sparse_stats(sp::SparseMatrixCSC)

Print some statistics for `sp`. The line
"number of non-zeros" counts structural non-zeros.
"""
function sparse_stats(sp::SparseMatrixCSC)
    padding = 18
    for (description, statistic) in (
        ("size", size(sp)),
        ("num. elements", length(sp)),
        ("num. non-zeros", nnz(sp)),
        ("sparsity", Printf.@sprintf("%.4e", sparsity(sp))),
        ("num. stored zeros", count(iszero, SparseArrays.nonzeros(sp))))
        println(rpad(description, padding), ": ", statistic)
    end
    return nothing
end

end # module SparseUtils

#  LocalWords:  numconnections countmap sp
