# Trying to use AbstractSparseMatrix now
# NO: At present, I use a Union type rather than subtyping (in abstract.jl), so that I don't have to
# implement the entire AbstractArray interface.
# const AbstractSparseMatrixUtils = Union{SparseMatrixCSC, SparseMatrixCOO}

export materialize

## SparseArrays uses `copy` for this. Using copy is confusing to me.
## Not very discoverable.
"""
    materialize(m::LinearAlgebra.Transpose)

`Transpose` is a thin wrapper around `m.parent`.
This function removes the wrapper and performs the
transpose. The documented way to do this is `copy(m)`.

Perform the material transform of the matrix `m.parent`.
This could call `copy` on `m`, or `permutedims` or
a specific tranpose routine on the parent. It seems
that `copy` is optimized, `permutedims` in sparsematrix.jl
falls back to generic methods (currently).
"""
function materialize(m::LinearAlgebra.Transpose)
    return copy(m)
#    return permutedims(m.parent)  not optimized for CSC
end
materialize(m) = m

"""
    density(sp::AbstractSparseMatrix)

Compute the density of `sp`, the fraction of elements that
are non-zero This assumes that there
are no stored values equal to zero.
"""
function density(sp::AbstractSparseMatrix)
    return SparseArrays.nnz(sp) / length(sp)
end

"""
    summarystats(sp::SparseMatrixCSC)

Print some statistics for `sp`. The line
"number of non-zeros" counts structural non-zeros.
"""
function summarystats(sp::AbstractSparseMatrix)
    padding = 18
    for (description, statistic) in (
        ("size", size(sp)),
        ("num. elements", length(sp)),
        ("num. non-zeros", nnz(sp)),
        ("density", Printf.@sprintf("%.4e", density(sp))),
        ("num. stored zeros", count(iszero, SparseArrays.nonzeros(sp))))
        println(rpad(description, padding), ": ", statistic)
    end
    return nothing
end
