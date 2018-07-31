using BenchmarkTools
using SparseArrays
using Revise; using SparseUtils

function genmats(m, n, density)
    spcsc = sprand(m, n, density)
    coo = SparseMatrixCOO(spcsc)
    return(spcsc, coo)
end

function buildmat(coo)
    a = zeros(eltype(coo.nzval), coo.m, coo.n)
    for j in 1:coo.n
        for i in 1:coo.m
            @show (i,j)
            a[i, j] = coo[i, j]
        end
    end
    a
end

slowArray(coo::SparseMatrixCOO) = reshape(collect(coo[i,j] for j in 1:ncols(coo) for i in 1:nrows(coo)), (nrows(coo), ncols(coo)))

nothing
