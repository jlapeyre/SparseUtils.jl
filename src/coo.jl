import SparseArrays:
    SparseMatrixCSC, AbstractSparseMatrix, nnz, nonzeros, nzvalview,
    dropzeros!, dropzeros

export
    SparseMatrixCOO

struct SparseMatrixCOO{Tv, Ti<:Integer} <: SparseArrays.AbstractSparseMatrix{Tv, Ti}
    m::Int
    n::Int
    Ir::Vector{Ti}
    Ic::Vector{Ti}
    nzval::Vector{Tv}
    function SparseMatrixCOO(m, n, Ir::Vector{Ti}, Ic::Vector{Ti}, nzval::Vector{Tv}) where {Ti <: Integer, Tv}
        length(Ir) == length(Ic) == length(nzval) || throw(DimensionMismatch("Ir, Ic, and nzval must have the same length."))
        (!isempty(Ir) && m < maximum(Ir)) &&
            throw(DimensionMismatch("Maximum row index $(maximum(Ir)) greater than row dimension $m."))
        (!isempty(Ic) && n < maximum(Ic)) && throw(DimensionMismatch("Maximum column index greater than column dimension."))
        new{Tv, Ti}(m, n, Ir, Ic, nzval)
    end
end

# Yes, we need this. It would have prevented a bug.
_splat_fields(S::SparseMatrixCOO) = (S.m, S.n, S.Ir, S.Ic, S.nzval)

Base.length(coo::SparseMatrixCOO) = prod(size(coo))
Base.size(coo::SparseMatrixCOO) = (coo.m, coo.n)
Base.eltype(coo::SparseMatrixCOO{Tv, Ti}) where {Tv, Ti} = Tv
Base.getindex(S::SparseMatrixCOO, i::Integer, j::Integer) = IJV.getindex(S.Ir, S.Ic, S.nzval, i, j)
Base.setindex!(S::SparseMatrixCOO, val, i::Integer, j::Integer) = IJV.setindex!(S.Ir, S.Ic, S.nzval, val, i, j)
SparseArrays.dropstored!(S::SparseMatrixCOO, i::Integer, j::Integer) = (IJV.dropstored!(S.Ir, S.Ic, S.nzval, i, j); S)
Base.issorted(S::SparseMatrixCOO) = IJV.issorted(S.Ir, S.Ic)

spzeros(args...; sparsetype=Type{SparseMatrixCOO}) = SparseMatrixCOO(IJV.spzeros(args...)...)
nnz(coo::SparseMatrixCOO) = length(coo.Ir)
nonzeros(coo::SparseMatrixCOO) = coo.nzval
ncols(coo::SparseMatrixCOO) = coo.n
nrows(coo::SparseMatrixCOO) = coo.m

# The @inline decorations and the hardcoding of functions in reduction increase
# speed dramatically. But, this may depend on Julia version.

@inline nzvalview(S::SparseMatrixCOO) = view(S.nzval, 1:nnz(S))
@inline Base.count(pred, S::SparseMatrixCOO) = NZVal.count(pred, S.m, S.n, nzvalview(S), length(S.nzval))

## FIXME: this is probably never called because SparseMatrixCOO is not a AbstractArray
Base._mapreduce(f, op, ::Base.IndexCartesian, S::SparseMatrixCOO) =
    NZVal._mapreduce(f, op, Base.IndexCartesian(), S.m, S.n, nzvalview(S), nnz(S))

@inline Base._mapreduce(f, op, S::SparseMatrixCOO) =
    NZVal._mapreduce(f, op, S.m, S.n, nzvalview(S), nnz(S))

for (fname, op) in ((:sum, :add_sum), (:maximum, :max), (:minimum, :min), (:prod, :mul_prod))
    @eval @inline (Base.$fname)(f::Callable, S::SparseMatrixCOO) = Base._mapreduce(f, Base.$op, S)
    @eval @inline (Base.$fname)(S::SparseMatrixCOO) = Base._mapreduce(identity, Base.$op, S)
end

@inline Base._mapreduce(f, op::Union{typeof(*), typeof(Base.mul_prod)}, S::SparseMatrixCOO{T}) where T =
    NZVal._mapreduce(f, op, S.m, S.n, nzvalview(S), nnz(S))

@inline Base.mapreduce(f, op, S::SparseMatrixCOO) = # Base._mapreduce(f, op, Base.IndexCartesian(), S)
    Base._mapreduce(f, op, S)

Base.copy(S::SparseMatrixCOO) =
    SparseMatrixCOO(S.m, S.n, Base.copy(S.Ir), Base.copy(S.Ic), Base.copy(S.nzval))

for fname in (:dropzeros, :dropzeros!)
    @eval ($fname)(S::SparseMatrixCOO) = (IJV.$fname)(S.Ir, S.Ic, S.nzval)
end

Base.:(==)(S::SparseMatrixCOO, V::SparseMatrixCOO) =
    S.m == V.m &&  S.n == V.n && S.Ir == V.Ir && S.nzval == V.nzval

# The following is somewhat inefficient. Iterating over nzvalview(coo) directly is better.
# But, defining `iterate` allows methods to be called with no more work.
Base.iterate(coo::SparseMatrixCOO, args...) = Base.iterate(SparseArrays.nzvalview(coo), args...)

### Conversion

function SparseMatrixCOO(spcsc::SparseMatrixCSC)
    SparseMatrixCOO(nrows(spcsc), ncols(spcsc), SparseArrays.findnz(spcsc)...)
end

SparseMatrixCSC(coo::SparseMatrixCOO) = SparseArrays.sparse(coo.Ir, coo.Ic, coo.nzval, coo.m, coo.n)
Base.Array(S::SparseMatrixCOO) = IJV.Array(_splat_fields(S)...)

prunecols(S::SparseMatrixCOO, min_entries; renumber=true) =
    SparseMatrixCOO(IJV.prunecols(_splat_fields(S)..., min_entries; renumber=renumber)...)
prunerows(S::SparseMatrixCOO, min_entries; renumber=true) =
    SparseMatrixCOO(IJV.prunerows(_splat_fields(S)..., min_entries; renumber=renumber)...)

## renumbercols, renumberrows
## FIXME: Factor out the builtin-only core.
# Wow, this is not super-readable
for (dim, Ind) in ((:cols, :Ic), (:rows, :Ir))
    renumberdim = Symbol(:renumber, dim)
    renumberdim! = Symbol(renumberdim, :!)
    _renumberdim! = Symbol(:_, renumberdim!)
    if dim == :cols
        call1 = :(SparseMatrixCOO(coo.m, last(newInd), coo.Ir, newInd, coo.nzval))
        call2 = :(SparseMatrixCOO(coo.m, last(coo.Ic), coo.Ir, coo.Ic, coo.nzval))
    else
        call1 = :(SparseMatrixCOO(maximum(newInd), coo.n, newInd, coo.Ic, coo.nzval))
        call2 = :(SparseMatrixCOO(maximum(coo.Ir), coo.n, coo.Ir, coo.Ic, coo.nzval))
    end
    @eval begin
        function $(renumberdim)(coo::SparseMatrixCOO{Tv,Ti}) where {Tv,Ti}
            newInd = Vector{Ti}(undef, length(coo.$Ind))
            $(_renumberdim!)(newInd, coo)
            return $call1
        end
        function $(renumberdim!)(coo::SparseMatrixCOO{Tv,Ti}) where {Tv,Ti}
            $(_renumberdim!)(coo.Ic, coo)
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

@doc """
    renumberrows(S::SparseMatrixCOO)

Remove rows that have no non-zero values
and renumber the rows.
""" renumberrows

@doc """
    renumberrows!(S::SparseMatrixCOO)

Remove rows that have no non-zero values
and renumber the rows.
""" renumberrows!

@doc """
    renumbercols(S::SparseMatrixCOO)

Remove cols that have no non-zero values
and renumber the cols.
""" renumbercols

@doc """
    renumbercols!(S::SparseMatrixCOO)

Remove cols that have no non-zero values
and renumber the cols.
""" renumbercols!

"""
    renumberrowscols(S::SparseMatrixCOO)

Equivalent to applying `renumberrow!` and `renumbercols`
"""
renumberrowscols(S::SparseMatrixCOO) = S |> renumberrows |> renumbercols

"""
    renumberrowscols!(S::SparseMatrixCOO)

Equivalent to applying `renumberrows!` and `renumbercols!`
"""
renumberrowscols!(S::SparseMatrixCOO) = S |> renumberrows! |> renumbercols!

for f in (:rotr90, :rotl90, :rot180)
    ijvf = Symbol(f, "!")
    @eval (Base.$f)(S::SparseMatrixCOO) = SparseMatrixCOO((IJV.$ijvf)(S.m, S.n, copy(S.Ir), copy(S.Ic), copy(S.nzval))...)
end

@doc """
    rotr90(S::SparseMatrixCOO)
rotate `S` clockwise 90 degrees.
""" rotr90

@doc """
    rotl90(S::SparseMatrixCOO)
rotate `S` counter-clockwise 90 degrees.
""" rotl90

@doc """
    rot180(S::SparseMatrixCOO)
rotate `S` 180 degrees.
""" rot180

"""
    permutedims!(S::SparseMatrixCOO)

Transpose `S` in place.
"""
Base.permutedims!(S::SparseMatrixCOO) = SparseMatrixCOO(IJV.permutedims!(_splat_fields(S)...)...)
#Base.permutedims!(S::SparseMatrixCOO) = SparseMatrixCOO(IJV.permutedims!(S.m, S.n, S.Ir, S.Ic, S.nzval)...)

"""
    permutedims(S::SparseMatrixCOO)

Return a transposed copy of `S`..
"""
Base.permutedims(S::SparseMatrixCOO) = permutedims!(copy(S))

# const _transpose_type = LinearAlgebra.Transpose{T, SparseMatrixCOO{T,V}} where V where T

function Base.copy(S::LinearAlgebra.Transpose{T, SparseMatrixCOO{T,V}}) where V where T
    return Base.permutedims(S.parent)
end

# We are now testing SparseMatrixCOO <: AbstractSparseMatrix
# The would require reproducing a lot of code that assumes the sparse matrix is an AbstractArray
#Base.transpose(S::SparseMatrixCOO) = throw(DomainError("Tranpose not implemented because `SparseMatrixCOO` is not an `AbstractArray`"))

function SparseArrays.sprand(args...; sparsetype=Type{SparseMatrixCOO})
    if sparsetype == SparseMatrixCOO
        return SparseMatrixCOO(IJV.sprand(args...)...)
    elseif sparsetype == SparseMatrixCSC
        return SparseArrays.sprand(args...)
    else
        throw(ArgumentError("Unsupported sparse matrix type"))
    end
end
