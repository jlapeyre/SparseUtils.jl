"""
    density(sp::AbstractSparseMatrix)

Compute the density of `sp`, the fraction of elements that
are non-zero. This assumes that there are no stored values equal to zero.
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
