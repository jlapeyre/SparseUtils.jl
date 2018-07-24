# SparseUtils

*Utilities for sparse matrices*

The source repository is [https://github.com/jlapeyre/SparseUtils.jl](https://github.com/jlapeyre/SparseUtils.jl).

## Contents

```@contents
```

## Index

```@index
```

## Functions

```@docs
sparsity
sparse_stats
nnz(sp::SparseMatrixCSC, colnum::Integer)
nnzcounts
c_to_julia_index!
c_to_julia_index
transpose(M::SparseMatrixCSC, lazy=true)
```