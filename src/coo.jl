export SparseMatrixCOO

struct SparseMatrixCOO{T <: Integer, W}
    I::Vector{T}
    J::Vector{T}
    V::Vector{W}
end

SparseMatrixCOO(spcsc::SparseArrays.SparseMatrixCSC) = SparseMatrixCOO(SparseArrays.findnz(spcsc)...)
SparseArrays.SparseMatrixCSC(coo::SparseMatrixCOO) = SparseArrays.sparse(coo.I, coo.J, coo.V)
