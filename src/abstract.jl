# At present, I use a Union type rather than subtyping (in abstract.jl), so that I don't have to
# implement the entire AbstractArray interface.
const AbstractSparseMatrixUtils = Union{SparseMatrixCSC, SparseMatrixCOO}

export materialize

"""
    materialize(m::LinearAlgebra.Transpose)

Perform the material transform of the matrix `m.parent`.

`Transpose` is a thin wrapper around `m.parent`.
This function removes the wrapper and performs the
transpose.
"""
function materialize(m::LinearAlgebra.Transpose)
    return material_transpose(m.parent)
end
materialize(m) = m

"""
    sparsity(sp::SparseMatrixCSC)

Compute the sparisty of `sp`. This assumes that there
are no stored values equal to zero.
"""
function sparsity(sp::AbstractSparseMatrixUtils)
    return SparseArrays.nnz(sp) / length(sp)
end

"""
    summarystats(sp::SparseMatrixCSC)

Print some statistics for `sp`. The line
"number of non-zeros" counts structural non-zeros.
"""
function summarystats(sp::AbstractSparseMatrixUtils)
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
