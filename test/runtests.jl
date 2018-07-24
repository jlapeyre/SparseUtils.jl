using SparseUtils
import SparseArrays
import SparseUtils.transpose
using Serialization
using Test

@testset "SparseUtils" begin
    sparse_array = open("sparse_array.dat", "r") do io
        deserialize(io)
    end
    expected_sparsity = 0.008095238095238095
    @test sparsity(sparse_array) == expected_sparsity

    expected_size = (30, 70)
    @test size(sparse_array) == expected_size
    @test transpose(sparse_array, lazy=false) |> size == (expected_size[2], expected_size[1])
    @test transpose(transpose(sparse_array, lazy=false), lazy=false) == sparse_array

    # nnz defined for columns
    @test sum(map(i -> SparseArrays.nnz(sparse_array, i),  1:size(sparse_array)[2])) == SparseArrays.nnz(sparse_array)
end
