var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "SparseUtils",
    "title": "SparseUtils",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SparseUtils-1",
    "page": "SparseUtils",
    "title": "SparseUtils",
    "category": "section",
    "text": "Utilities for sparse matricesThe source repository is https://github.com/jlapeyre/SparseUtils.jl."
},

{
    "location": "index.html#Contents-1",
    "page": "SparseUtils",
    "title": "Contents",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Index-1",
    "page": "SparseUtils",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#SparseUtils.sparsity",
    "page": "SparseUtils",
    "title": "SparseUtils.sparsity",
    "category": "function",
    "text": "sparsity(sp::SparseMatrixCSC)\n\nCompute the sparisty of sp. This assumes that there are no stored values equal to zero.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.sparse_stats",
    "page": "SparseUtils",
    "title": "SparseUtils.sparse_stats",
    "category": "function",
    "text": "sparse_stats(sp::SparseMatrixCSC)\n\nPrint some characteristics of sp.\n\n\n\n\n\n"
},

{
    "location": "index.html#Base.nnz",
    "page": "SparseUtils",
    "title": "Base.nnz",
    "category": "function",
    "text": "nnz(sp::SparseMatrixCSC, colnum::Integer)\n\nReturn the number of structurally non-zero elements in the colnumth column of sp.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.c_to_julia_index!",
    "page": "SparseUtils",
    "title": "SparseUtils.c_to_julia_index!",
    "category": "function",
    "text": "c_to_julia_index!(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC\n\nConvert a sparse matrix with zero-based indices to a SparseMatrixCSC. colptr0 and rowval0 are altered in place. nzval is not copied.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.c_to_julia_index",
    "page": "SparseUtils",
    "title": "SparseUtils.c_to_julia_index",
    "category": "function",
    "text": "c_to_julia_index(colptr0, rowval0, nzval)::SparseArrays.SparseMatrixCSC\n\nConvert a sparse matrix with zero-based indices to a SparseMatrixCSC. colptr0 and rowval0 are copied and altered. nzval is not copied.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.transpose_concrete",
    "page": "SparseUtils",
    "title": "SparseUtils.transpose_concrete",
    "category": "function",
    "text": "transpose_concrete(M::SparseMatrixCSC)\n\nReturn the concrete transpose of the sparse matrix M. That is, return a new sparse matrix, ratherM wrapped in Transpose.\n\n\n\n\n\n"
},

{
    "location": "index.html#Functions-1",
    "page": "SparseUtils",
    "title": "Functions",
    "category": "section",
    "text": "sparsity\nsparse_stats\nnnz\nc_to_julia_index!\nc_to_julia_index\ntranspose_concrete"
},

]}
