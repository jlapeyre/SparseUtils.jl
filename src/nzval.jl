module NZVal

## There is some code in sparsematrix.jl that depends only on the
## number of rows and columns, the collection of nonzero values,
## and possibly the number of nonzero values.
## In particular, the code does not depend on the column and row
## information in colptr and rowval.
## This code is factored out here. It can be used directly in
## CSC, CSR, COO, and possibly other sparse matrix representations.

## _mapreducezeros in sparsematrix.jl does not depend on any sparse type.
## So, we do not need to reproduce it here.

import SparseArrays: _mapreducezeros

# for sym in (:length, :count, :mapreduce)
#     Core.eval(SparseArrays, :(function $sym end))
# end

# import SparseArrays: length, count, mapreduce

mnlength(m, n) = m * n

# count(pred, S::SparseMatrixCSC) = count(pred, nzvalview(S)) + pred(zero(eltype(S)))*(prod(size(S)) - nnz(S))
_count(pred, m, n, nzval::AbstractArray{T}, numnz=length(nzval)) where T =
    Base.count(pred, nzval) + pred(zero(T)) * (prod(m, n) - numnz)

#function Base._mapreduce(f, op, ::Base.IndexCartesian, A::SparseMatrixCSC{T}) where T
function _mapreduce(f, op, ::Base.IndexCartesian, m, n, nzval::AbstractArray{T}, numnz) where T
    return _mapreduce(f, op, n, n, nzval, numnz)
end

function _mapreduce(f, op, m, n, nzval::AbstractArray{T}, numnz) where T
    z = numnz
    nlength = mnlength(m, n)
    if z == 0
        if nlength == 0
            Base.mapreduce_empty(f, op, T)
        else
            _mapreducezeros(f, op, T, nlength-z-1, f(zero(T)))
        end
    else
        _mapreducezeros(f, op, T, nlength-z, Base._mapreduce(f, op, nzval)) # was nzvalview(A)
    end
end

# nzvalview(A)) --> nzval
# function Base._mapreduce(f, op::typeof(*), A::SparseMatrixCSC{T}) where T
function _mapreduce(f, op::Union{typeof(*), typeof(Base.mul_prod)},
                      m, n, nzval::AbstractArray, numnz=length(nzval))
    nzeros = mnlength(m, n) - numnz
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
