using SparseUtils
import SparseArrays
import SparseArrays: SparseMatrixCSC
import SparseUtils: transpose, materialize
import SparseUtils:  SparseMatrixCOO
using Serialization
using Test

let # typeof(sparse_array) = SparseMatrixCSC
    sparse_array = open("sparse_array.dat", "r") do io
        deserialize(io)
    end
    expected_size = (30, 70)
    expected_sum = 1.6277046559059591
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

    let sparse_array_coo = SparseMatrixCOO(sparse_array)

        @testset "conversion" begin
            @test SparseMatrixCSC(sparse_array_coo) == sparse_array
        end

        @testset "sum" begin
            @test sum(sparse_array) == sum(sparse_array_coo)
        end

        @testset "prod" begin
            @test prod(sparse_array) == prod(sparse_array_coo)
            zerotoone(x) = x == 0 ? one(x) : x
            @test prod(zerotoone, sparse_array) == prod(zerotoone, sparse_array_coo)
        end
    end

#    @test isapprox(sum(sparse_array), expected_sum)
end
