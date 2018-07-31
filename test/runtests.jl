using SparseUtils
import SparseArrays
import SparseArrays: SparseMatrixCSC
import SparseUtils: materialize
import SparseUtils: SparseMatrixCOO
import LinearAlgebra
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
        @test length(sparse_array) == prod(size(sparse_array))
    end

    @testset "density" begin
        expected_density = 0.008095238095238095
        @test density(sparse_array) == expected_density
    end

    @testset "transpose" begin
        @test transpose(sparse_array) |> materialize |> size == (expected_size[2], expected_size[1])
        @test sparse_array |> transpose |> materialize |> transpose |> materialize == sparse_array
    end

    @testset "nnz" begin
        # nnz defined for columns
        @test sum(map(i -> SparseArrays.nnz(sparse_array, i),  1:size(sparse_array)[2])) == SparseArrays.nnz(sparse_array)
    end

    @testset "construction" begin
        S = SparseArrays.sparse([1], [1], [1], 1, 1; sparsetype=SparseMatrixCOO)
        @test isa(S, SparseMatrixCOO{Int,Int})
        @test size(S) == (1, 1)
        @test length(S) == 1
        S1 = SparseArrays.sparse([5, 7], [2, 1], [1.0, 2.0], 10, 10; sparsetype=SparseMatrixCOO)
        @test size(S1) == (10, 10)
    end

    let sparse_array_coo = SparseMatrixCOO(sparse_array)

        @testset "conversion" begin
            @test SparseMatrixCSC(sparse_array_coo) == sparse_array
        end

        @testset "transpose" begin
            pdcoo = permutedims(sparse_array_coo)
            @test copy(transpose(sparse_array_coo)) == pdcoo
            @test isa(transpose(sparse_array_coo), LinearAlgebra.Transpose)
            @test permutedims(pdcoo) == sparse_array_coo
            @test sparse_array_coo |> transpose |> transpose == sparse_array_coo
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
