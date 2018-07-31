# SparseUtils

*Utilities for sparse matrices*

The source repository is [https://github.com/jlapeyre/SparseUtils.jl](https://github.com/jlapeyre/SparseUtils.jl).

This package implements the type `SparseMatrixCOO`, several functions for type `SparseMatrixCSC`,
and a few functions for `AbstractSparseMatrix`. Submodules `NZVal` and `IJV` contain the generic
algorithms, which are wrapped by the sparse matrix types.

Several functions are not yet documented.

## Contents

```@contents
```

## Index

```@index
```

## SparseMatrixCOO

```@docs
SparseUtils.SparseMatrixCOO
SparseUtils.sparse
SparseUtils.issorted
SparseUtils.spzeros
SparseUtils.prunecols
SparseUtils.prunerows
```

## SparseMatrixCSC

```@docs
SparseUtils.nnz(sp::SparseMatrixCSC, colnum::Integer)
SparseUtils.nnzcounts
SparseUtils.c_to_julia_index!
SparseUtils.c_to_julia_index
```

<!-- SparseUtils.transpose -->

## AbstractSparseMatrix

```@docs
SparseUtils.density
SparseUtils.summarystats
```

## NZVal

```@docs
SparseUtils.NZVal
```

## IJV

```@docs
SparseUtils.IJV
SparseUtils.IJV.issorted
SparseUtils.IJV.findindex
SparseUtils.IJV.getindex
```
