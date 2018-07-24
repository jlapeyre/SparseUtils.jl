module SparseUtils

export c_to_julia_index, c_to_julia_index!, sparsity,
    nnzcounts, sparse_stats, nrows, ncols,
    prunecols!, dropzerocols

import SparseArrays
import SparseArrays: SparseMatrixCSC, nnz
import DataUtils
import StatsBase
import Printf

"""
    c_to_julia_index!(colptr, rowval, nzval)::SparseArrays.SparseMatrixCSC

Convert a sparse matrix with zero-based indices to a `SparseMatrixCSC`.
`colptr` and `rowval` are altered in place. `nzval` is neither copied
nor altered.
"""
function c_to_julia_index!(colptr, rowval, nzval)::SparseArrays.SparseMatrixCSC
    colptr .= colptr .+ 1
    rowval .= rowval .+ 1
    n = length(colptr)-1
    m = maximum(rowval)
    return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end

"""
    c_to_julia_index(colptr, rowval, nzval)::SparseArrays.SparseMatrixCSC

Convert a sparse matrix with zero-based indices to a `SparseMatrixCSC`.
`colptr` and `rowval` are not altered, but rather copied. `nzval` is neither copied, nor altered.
"""
function c_to_julia_index(colptr, rowval, nzval)::SparseArrays.SparseMatrixCSC
    c_to_julia_index!(deepcopy(colptr), deepcopy(rowval), nzval)
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
    transpose(M::SparseMatrixCSC; lazy=true)

Return the lazy transpose of `M` if `lazy=true` and the material
transpose, of type `SparseMatrixCSC`, otherwise.
"""
function Base.transpose(M::SparseMatrixCSC; lazy=true)
    return lazy ? SparseArrays.Transpose(M) : SparseArrays.ftranspose(M, identity)
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
    nnzcounts(sp::SparseMatrixCSC; rev=true, byvalue=true)

Return a sorted countmap of the number of structural non-zeros
in the columns of `sp`. The keys controlling the sorting are
passed to `StatsBase.sort`.

# Examples
In this example, there are 132781 columns with one non-zero element,
32096 columns with two, etc.
```julia
julia> conns = nnzcounts(sp);

julia> collect(Iterators.take(conns,3))
3-element Array{Pair{Int64,Int64},1}:
 1 => 132781
 2 => 32096
 3 => 17478
```
"""
function nnzcounts(sp::SparseMatrixCSC; rev=true, byvalue=true)
    spcountmap = DataUtils.countmap(nnzcols(sp), datatype=Int)
    sorted_spcountmap = StatsBase.sort(spcountmap, rev=rev, byvalue=rev)
    return sorted_spcountmap
end

# These may have been in SparseArrays, but were removed.
# They are very convenient at the command line.
nrows(sp::SparseMatrixCSC) = size(sp)[1]
ncols(sp::SparseMatrixCSC) = size(sp)[2]

"""
    nnzcols(sp::SparseMatrixCSC)

Return an iterator over the number of non-zero elements
in each column of `sp`.
"""
function nnzcols(sp::SparseMatrixCSC)
    (nnz(sp, i) for i in 1:ncols(sp))
end

function prunecols!(sp::SparseMatrixCSC, min_connections)
    connection_flags = [nnz(sp, j) >= min_connections for j in 1:ncols(sp)]
    SparseArrays.fkeep!(sp, (i, j, v) -> connection_flags[j], true)
    return sp
end

function dropzerocols(sp::SparseMatrixCSC)
    newcolptr = unique(sp.colptr)
    n = length(newcolptr) - 1
    SparseMatrixCSC(sp.m, n, newcolptr, sp.rowval, sp.nzval)
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

# function prune(sp::SparseMatrixCSC)
#     (I, J, V) = SparseArrays.findnz(sp)
# end

end # module SparseUtils

#  LocalWords:  nnzcounts countmap sp
