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
