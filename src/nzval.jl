"""
    module NZVal

For sparse matrices, algorithms that depend only on nzval, m, and n.

This code is generic, for types in `Core`.

The symbols in `NZVal` should *not* be exported, not even for
convenience. They conflict with functions in `Base`.
"""
module NZVal

## This module contains code that depends only on a container
## of non-zero values and the dimensions of the matrix. The
## code references only builtin types.
##
## There is some code in sparsematrix.jl that depends only on the
## number of rows and columns, the collection of nonzero values,
## and possibly the number of nonzero values.
## In particular, the code does not depend on the column and row
## information in colptr and rowval.
## This code is factored out here. It can be used directly in
## CSC, CSR, COO, and possibly other sparse matrix representations.

##  _mapreducezeros in sparsematrix.jl does not depend on any sparse type.
##  So, we import `_mapreducezeros` rather than reproducing it here.
import SparseArrays: _mapreducezeros

@inline numzeros(m, n, numnz) = m * n - numnz

count(pred, m::Integer, n::Integer, nzval::AbstractArray{T}, numnz::Integer) where T =
      Base.count(pred, nzval) +  pred(zero(T)) * numzeros(m, n, numnz)

# sparsematrix.jl: function Base._mapreduce(f, op, ::Base.IndexCartesian, A::SparseMatrixCSC{T}) where T
function _mapreduce(f, op, ::Base.IndexCartesian, m, n, nzval::AbstractArray, numnz)
    return _mapreduce(f, op, n, n, nzval, numnz)
end

function _mapreduce(f, op, m, n, nzval::AbstractArray{T}, numnz) where T
    z = numnz
    nzeros = numzeros(m, n, numnz)
    nlength = m * n
    if z == 0
        if nlength == 0
            Base.mapreduce_empty(f, op, T)
        else
            _mapreducezeros(f, op, T, nzeros - 1, f(zero(T)))
        end
    else
        _mapreducezeros(f, op, T, nzeros, Base._mapreduce(f, op, nzval)) # was nzvalview(A)
    end
end

# sparsematrix.jl: function Base._mapreduce(f, op::typeof(*), A::SparseMatrixCSC{T}) where T
function _mapreduce(f, op::Union{typeof(*), typeof(Base.mul_prod)},
                      m, n, nzval::AbstractArray, numnz=length(nzval))
    nzeros = numzeros(m, n, numnz)
    if nzeros == 0
        # No zeros, so don't compute f(0) since it might throw
        Base._mapreduce(f, op, nzval)
    else
        T = eltype(nzval)
        v = f(zero(T))^(nzeros)
        # Bail out early if initial reduction value is zero
        v == zero(T) ? v : v * Base._mapreduce(f, op, nzval)
    end
end

end # module
