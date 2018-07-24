# SparseUtils

*Utilities for sparse matrices*

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://jlapeyre.github.io/SparseUtils.jl/latest)
Linux, OSX: [![Build Status](https://travis-ci.org/jlapeyre/SparseUtils.jl.svg?branch=master)](https://travis-ci.org/jlapeyre/SparseUtils.jl)
&nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/SparseUtils.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/sparseutils-jl)
&nbsp; &nbsp; &nbsp;
[![Coverage Status](https://coveralls.io/repos/jlapeyre/SparseUtils.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jlapeyre/SparseUtils.jl?branch=master)
[![codecov.io](http://codecov.io/github/jlapeyre/SparseUtils.jl/coverage.svg?branch=master)](http://codecov.io/github/jlapeyre/SparseUtils.jl?branch=master)

See the documentation pages or the doc strings for 

* `nnz` extended to operate on columns
* `c_to_julia_index`, `c_to_julia_index!` import zero-based-indexing data
* `sparsity` compute sparsity
* `transpose_concrete` return another `SparseMatrixCSC`
* `sparse_stats` print some measurements
* `numconnections` count map over columns of number of non-zeros
