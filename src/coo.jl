import SparseArrays.SparseMatrixCSC
import SparseArrays.AbstractSparseMatrix

export
    SparseMatrixCOO,
    renumbercols!

# At present, I use a Union type rather than subtyping (in abstract.jl), so that I don't have to
# implement the entire AbstractArray interface.
struct SparseMatrixCOO{T <: Integer, W} # <: SparseArrays.AbstractSparseMatrix{T, W}
    m::Int
    n::Int
    I1::Vector{T}
    J::Vector{T}
    V::Vector{W}
    function SparseMatrixCOO(m, n, I1::Vector{T}, J::Vector{T}, V::Vector{W}) where {T <: Integer, W}
        length(I1) == length(J) == length(V) || throw(DimensionMismatch("I, J, and V must have the same length."))
        m < maximum(I1) && throw(DimensionMismatch("Maximum row index $(maximum(I1)) greater than row dimension $m."))
        n < maximum(J) && throw(DimensionMismatch("Maximum column index greater than column dimension."))
        new{T,W}(m, n, I1, J, V)
    end
end

SparseArrays.nnz(coo::SparseMatrixCOO) = length(coo.I1)
SparseArrays.nonzeros(coo::SparseMatrixCOO) = coo.V
ncols(coo::SparseMatrixCOO) = coo.n
nrows(coo::SparseMatrixCOO) = coo.m
Base.length(coo::SparseMatrixCOO) = coo.n * coo.m
Base.size(coo::SparseMatrixCOO) = (coo.m, coo.m)
# Base.getindex(coo::SparseMatrixCOO, i::Integer, j::Integer) =

### Conversion

function SparseMatrixCOO(spcsc::SparseMatrixCSC)
    SparseMatrixCOO(nrows(spcsc), ncols(spcsc), SparseArrays.findnz(spcsc)...)
end

function SparseMatrixCSC(coo::SparseMatrixCOO)
    return SparseArrays.sparse(coo.I1, coo.J, coo.V, coo.m, coo.n)
end

# Wow, this is not super-readable
for (dim, Ind) in ((:cols, :J), (:rows, :I1))
    renumberdim = Symbol(:renumber, dim)
    renumberdim! = Symbol(renumberdim, :!)
    _renumberdim! = Symbol(:_, renumberdim!)
    if dim == :cols
        call1 = :(SparseMatrixCOO(coo.m, last(newInd), coo.I1, newInd, coo.V))
        call2 = :(SparseMatrixCOO(coo.m, last(coo.J), coo.I1, coo.J, coo.V))
    else
        call1 = :(SparseMatrixCOO(last(newInd), coo.n, newInd, coo.J, coo.V))
        call2 = :(SparseMatrixCOO(last(coo.I1), coo.n, coo.I1, coo.J, coo.V))
    end
    @eval begin
        function $(renumberdim)(coo::SparseMatrixCOO{T,W}) where {T,W}
            newInd = Vector{T}(undef, length(coo.$Ind))
            $(_renumberdim!)(newInd, coo)
            return $call1
        end
        function $(renumberdim!)(coo::SparseMatrixCOO{T,W}) where {T,W}
            $(_renumberdim!)(coo.J, coo)
            return $call2
        end
        function $(_renumberdim!)(newInd, coo::SparseMatrixCOO{T,W}) where {T,W}
            oldIndlist = sort!(unique(coo.$Ind))
            for j in 1:length(coo.$Ind)
                newInd[j] = searchsortedfirst(oldIndlist, coo.$Ind[j])
            end
            return nothing
        end
    end
end
