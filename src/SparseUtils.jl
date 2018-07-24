module SparseUtils

export c_to_julia_index, c_to_julia_index!, sparsity,
       transpose_concrete, sparse_stats

import SparseArrays
import SparseArrays: SparseMatrixCSC, nnz

"""
    c_to_julia_index!(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC

Convert a sparse matrix with zero-based indices to a `SparseMatrixCSC`.
`colptr0` and `rowval0` are altered in place. `nzval` is not copied.
"""
function c_to_julia_index!(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC
    colptr0 .+= 1
    rowval0 .+ 1
    m = length(colptr0)-1
    n = maximum(rowval0)
    return SparseArrays.SparseMatrixCSC(n, m, colptr0, rowval0, nzval)
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
return a new sparse matrix, rather`M` wrapped in `Transpose`.
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
    sparse_stats(sp::SparseMatrixCSC)

Print some characteristics of `sp`.
"""
function sparse_stats(sp::SparseMatrixCSC)
    for s in (
        string("sparse matrix size: size(sp) = ", size(sp)),
        string("number of elements: length(sp) = ", length(sp)),
        string("number of non-zero elements: nnz(sp) = ", SparseArrays.nnz(sp)),
        string("sparsity: ", sparsity(sp)))

        println(s)
    end
    return nothing
end

end # module SparseUtils
