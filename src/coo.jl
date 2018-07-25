struct SparseMatrixCOO{T,V} where {T<:Integer, W}
    I::T
    J::T
    V::W
end
