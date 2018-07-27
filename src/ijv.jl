module IJV

## Algorithms for coo representation of sparse matrices
## This module only makes use of builtin types.

function dropzeros!(I, J, V::AbstractArray{T}) where T
    mask = V .== zero(T)
    map(x -> deleteat!(x, mask), (I, J, V))
end
dropzeros(I, J, V::AbstractArray) = dropzeros!(copy(I), copy(J), copy(V))

# This is copied from sparsematrix.jl, modified only to change the order of the arguments
function rot180!(m, n, I, J, V)
    for i=1:length(I)
        I[i] = m - I[i] + 1
        J[i] = n - J[i] + 1
    end
    return (m, n, I, J, V)
end

function rotr90!(m, n, I, J, V)
    #old col inds are new row inds
    for i=1:length(I)
        I[i] = m - I[i] + 1
    end
    return (n, m, J, I, V)
end

function rotl90!(m, n, I, J, V)
    #old row inds are new col inds
    for i=1:length(J)
        J[i] = n - J[i] + 1
    end
    return (n, m, J, I, V)
end

end # module
