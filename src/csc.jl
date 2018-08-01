## The code in this file *does* depend on SparseMatrixCSC
## The algorithms could be factored out.

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

# We should use a doc like this.
# """
#     transpose(M::SparseMatrixCSC)

# Return the lazy transpose of `M`, which is a thin wrapper
# around `M`. To get the material
# transpose, which is another matrix of type `SparseMatrixCSC`,
# use `materialize(transpose(M))`.
# """

## FIXME: probably should remove this.
material_transpose(M::SparseMatrixCSC) = SparseArrays.ftranspose(M, identity)

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
    spcountmap = SparseUtils.countmap(nnzcols(sp), datatype=Int)
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
nnzcols(sp::SparseMatrixCSC) = (nnz(sp, i) for i in 1:ncols(sp))

"""
    prunecols!(sp::SparseMatrixCSC, min_connections)

For each column, set all elements of the column equal to zero if the number of non-zero
elements is less than `min_connections`. `sp` is modified in place. Here "zero" means
structural zeros.
"""
function prunecols!(sp::SparseMatrixCSC, min_connections)
    connection_flags = [nnz(sp, j) >= min_connections for j in 1:ncols(sp)] # This must not change during fkeep!
    SparseArrays.fkeep!(sp, (i, j, v) -> connection_flags[j], true)
    return sp
end

"""
    hasemptycols(sp::SparseMatrixCSC)

Return true if `sp` contains columns for which all elements are structural
zeros.

`prunecols!` will create such columns.
"""
hasemptycols(sp::SparseMatrixCSC) = length(sp.colptr) != length(unique(sp.colptr))

"""
    renumbercols(sp::SparseMatrixCSC)

If `sp` has "empty" columns as determined by `hasemptycols`, then remove these
columns and renumber the remaining columns. A new `SparseMatrixCSC` is returned.
Modifying elements of `.rowval` and `.nzval` in the output matrix will change
these values in `sp`.
"""
function renumbercols(sp::SparseMatrixCSC)
    (newcolptr, n) = _renumbercols(sp)
    SparseMatrixCSC(sp.m, n, newcolptr, sp.rowval, sp.nzval)
end

function _renumbercols(sp::SparseMatrixCSC)
    newcolptr = unique(sp.colptr)
    n = length(newcolptr) - 1
    return (newcolptr, n)
end

"""
    hasemptyrows(sp::SparseMatrixCSC)

Return true if `sp` contains rows for which all elements are structural
zeros.

`prunecols!` may create such rows.
"""
function hasemptyrows(sp::SparseMatrixCSC)
    rowlist = unique(sp.rowval)
    return length(rowlist) != maximum(rowlist)
end

"""
    renumberrows(sp::SparseMatrixCSC)

If `sp` has "empty" rows as determined by `hasemptycols`, then remove these
columns and renumber the remaining columns. A new `SparseMatrixCSC` is returned.
Modifying elements of `.rowval` and `.nzval` in the output matrix will change
these values in `sp`.
"""
function renumberrows(sp::SparseMatrixCSC)
    (renumbered_rowval, m) = _renumberrows(sp)
    SparseMatrixCSC(m, sp.n, sp.colptr, renumbered_rowval, sp.nzval)
end

function _renumberrows(sp::SparseMatrixCSC)
    sorted_rowval = sort!(unique(sp.rowval))
    rowval_type = eltype(sp.rowval)
    renumbered_rowval = rowval_type[searchsortedfirst(sorted_rowval, val) for val in sp.rowval]
    m = length(sorted_rowval)
    return (renumbered_rowval, m)
end

"""
    renumberrowscols(sp::SparseMatrixCSC)

Perform both `renumberrows` and `renumbercols` on `sp`
"""
function renumberrowscols(sp::SparseMatrixCSC)
    (newcolptr, n) = _renumbercols(sp)
    (renumbered_rowval, m) = _renumberrows(sp)
    return SparseMatrixCSC(m, n, newcolptr, renumbered_rowval, sp.nzval)
end
