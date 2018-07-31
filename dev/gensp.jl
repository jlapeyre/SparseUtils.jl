using SparseArrays
using SparseUtils
using Serialization

function create_serialize_sparse_array(nr=30, nc=70, density=0.01)
    sparse_array = SparseArrays.sprandn(nr, nc, density)
    open("sparse_array.dat", "w") do io
        serialize(io, sparse_array)
    end
    return nothing
end

function deserialize_sparse_array()
    return open("sparse_array.dat", "r") do io
        res = deserialize(io)
    end
end

