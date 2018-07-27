import SparseArrays:
    SparseMatrixCSC, AbstractSparseMatrix, nnz, nonzeros, nzvalview,
    dropzeros!, dropzeros

export
    SparseMatrixCOO,
    renumbercols!

struct SparseMatrixCOO{Tv, Ti<:Integer} # <: SparseArrays.AbstractSparseMatrix{T, W}
    m::Int
    n::Int
    Ir::Vector{Ti}
    Ic::Vector{Ti}
    nzval::Vector{Tv}
    function SparseMatrixCOO(m, n, Ir::Vector{Ti}, Ic::Vector{Ti}, nzval::Vector{Tv}) where {Ti <: Integer, Tv}
        length(Ir) == length(Ic) == length(nzval) || throw(DimensionMismatch("Ir, Ic, and nzval must have the same length."))
        m < maximum(Ir) && throw(DimensionMismatch("Maximum row index $(maximum(Ir)) greater than row dimension $m."))
        n < maximum(Ic) && throw(DimensionMismatch("Maximum column index greater than column dimension."))
        new{Tv, Ti}(m, n, Ir, Ic, nzval)
    end
end

Base.eltype(coo::SparseMatrixCOO{Tv, Ti}) where {Tv, Ti} = Tv
nnz(coo::SparseMatrixCOO) = length(coo.Ir)
nonzeros(coo::SparseMatrixCOO) = coo.nzval
ncols(coo::SparseMatrixCOO) = coo.n
nrows(coo::SparseMatrixCOO) = coo.m
Base.length(coo::SparseMatrixCOO) = coo.n * coo.m
Base.size(coo::SparseMatrixCOO) = (coo.m, coo.m)

# The @inline decorations and the hardcoding of functions for reduction increase
# speed dramatically. But, this may depend on Julia version.

nzvalview(S::SparseMatrixCOO) = view(S.nzval, 1:nnz(S))
@inline Base.count(pred, S::SparseMatrixCOO) = NZVal._count(pred, S.m, S.n, nzvalview(S), length(S.nzval))

## FIXME: this is probably never called because SparseMatrixCOO is not a AbstractArray
Base._mapreduce(f, op, ::Base.IndexCartesian, S::SparseMatrixCOO) =
    NZVal._mapreduce(f, op, Base.IndexCartesian(), S.m, S.n, nzvalview(S), nnz(S))

@inline Base._mapreduce(f, op, S::SparseMatrixCOO) =
    NZVal._mapreduce(f, op, S.m, S.n, nzvalview(S), nnz(S))

for (fname, op) in ((:sum, :add_sum), (:maximum, :max), (:minimum, :min), (:prod, :mul_prod))
    @eval @inline (Base.$fname)(f, S::SparseMatrixCOO) = Base._mapreduce(f, Base.$op, S)
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

function SparseMatrixCSC(coo::SparseMatrixCOO)
    return SparseArrays.sparse(coo.Ir, coo.Ic, coo.nzval, coo.m, coo.n)
end

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
renumberrowscols(S::SparseMatrixCOO) = S |> renumberrows |> renumbercols

for f in (:rotr90, :rotl90, :rot180)
    ijvf = Symbol(f, "!")
    @eval (Base.$f)(S::SparseMatrixCOO) = SparseMatrixCOO((IJV.$ijvf)(S.m, S.m, copy(S.Ir), copy(S.Ic), copy(S.nzval))...)
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
