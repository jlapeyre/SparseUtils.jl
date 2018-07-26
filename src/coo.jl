import SparseArrays.SparseMatrixCSC
import SparseArrays.AbstractSparseMatrix

export
    SparseMatrixCOO,
    renumbercols!

# At present, I use a Union type rather than subtyping (in abstract.jl), so that I don't have to
# implement the entire AbstractArray interface.
struct SparseMatrixCOO{Tv, Ti<:Integer} # <: SparseArrays.AbstractSparseMatrix{T, W}
    m::Int
    n::Int
    I1::Vector{Ti}
    J::Vector{Ti}
    nzval::Vector{Tv}
    function SparseMatrixCOO(m, n, I1::Vector{Ti}, J::Vector{Ti}, nzval::Vector{Tv}) where {Ti <: Integer, Tv}
        length(I1) == length(J) == length(nzval) || throw(DimensionMismatch("I, J, and nzval must have the same length."))
        m < maximum(I1) && throw(DimensionMismatch("Maximum row index $(maximum(I1)) greater than row dimension $m."))
        n < maximum(J) && throw(DimensionMismatch("Maximum column index greater than column dimension."))
        new{Tv, Ti}(m, n, I1, J, nzval)
    end
end

Base.eltype(coo::SparseMatrixCOO{Tv, Ti}) where {Tv, Ti} = Tv
SparseArrays.nnz(coo::SparseMatrixCOO) = length(coo.I1)
SparseArrays.nonzeros(coo::SparseMatrixCOO) = coo.nzval
ncols(coo::SparseMatrixCOO) = coo.n
nrows(coo::SparseMatrixCOO) = coo.m
Base.length(coo::SparseMatrixCOO) = coo.n * coo.m
Base.size(coo::SparseMatrixCOO) = (coo.m, coo.m)

SparseArrays.nzvalview(S::SparseMatrixCOO) = view(S.nzval, 1:nnz(S))
Base.count(pred, S::SparseMatrixCOO) = Base.count(pred, SparseArrays.nzvalview(S)) +
    pred(zero(eltype(S)))*(prod(size(S)) - nnz(S))

# Several _mapreduce methods are needed. Many are exactly the same as for CSC,
# but are typed for SparseMatrixCSC. So, we can't call those.
# In the meantime, many methods, such as prod, silently return the wrong result.

# The following is somewhat inefficient. Iterating over nzvalview(coo) directly is better.
# But, defining iterate allows many methods to be called with no more work.
Base.iterate(coo::SparseMatrixCOO, args...) = Base.iterate(SparseArrays.nzvalview(coo), args...)

### Conversion

function SparseMatrixCOO(spcsc::SparseMatrixCSC)
    SparseMatrixCOO(nrows(spcsc), ncols(spcsc), SparseArrays.findnz(spcsc)...)
end

function SparseMatrixCSC(coo::SparseMatrixCOO)
    return SparseArrays.sparse(coo.I1, coo.J, coo.nzval, coo.m, coo.n)
end

# Wow, this is not super-readable
for (dim, Ind) in ((:cols, :J), (:rows, :I1))
    renumberdim = Symbol(:renumber, dim)
    renumberdim! = Symbol(renumberdim, :!)
    _renumberdim! = Symbol(:_, renumberdim!)
    if dim == :cols
        call1 = :(SparseMatrixCOO(coo.m, last(newInd), coo.I1, newInd, coo.nzval))
        call2 = :(SparseMatrixCOO(coo.m, last(coo.J), coo.I1, coo.J, coo.nzval))
    else
        call1 = :(SparseMatrixCOO(last(newInd), coo.n, newInd, coo.J, coo.nzval))
        call2 = :(SparseMatrixCOO(last(coo.I1), coo.n, coo.I1, coo.J, coo.nzval))
    end
    @eval begin
        function $(renumberdim)(coo::SparseMatrixCOO{Tv,Ti}) where {Tv,Ti}
            newInd = Vector{T}(undef, length(coo.$Ind))
            $(_renumberdim!)(newInd, coo)
            return $call1
        end
        function $(renumberdim!)(coo::SparseMatrixCOO{Tv,Ti}) where {Tv,Ti}
            $(_renumberdim!)(coo.J, coo)
            return $call2
        end
        function $(_renumberdim!)(newInd, coo::SparseMatrixCOO{Tv,Ti}) where {Tv,Ti}
            oldIndlist = sort!(unique(coo.$Ind))
            for j in 1:length(coo.$Ind)
                newInd[j] = searchsortedfirst(oldIndlist, coo.$Ind[j])
            end
            return nothing
        end
    end
end
