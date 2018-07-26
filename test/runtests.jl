using SparseUtils
import SparseArrays
import SparseUtils: transpose, materialize
using Serialization
using Test

let
    sparse_array = open("sparse_array.dat", "r") do io
        deserialize(io)
    end
    expected_size = (30, 70)

    @testset "measurements" begin
        @test size(sparse_array) == expected_size
    end

    @testset "sparsity" begin
        expected_sparsity = 0.008095238095238095
        @test sparsity(sparse_array) == expected_sparsity
    end

    @testset "transpose" begin
        @test transpose(sparse_array) |> materialize |> size == (expected_size[2], expected_size[1])
        @test sparse_array |> transpose |> materialize |> transpose |> materialize == sparse_array
    end

    @testset "nnz" begin
        # nnz defined for columns
        @test sum(map(i -> SparseArrays.nnz(sparse_array, i),  1:size(sparse_array)[2])) == SparseArrays.nnz(sparse_array)
    end
end
