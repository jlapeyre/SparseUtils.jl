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
    "text": "Utilities for sparse matricesThe source repository is https://github.com/jlapeyre/SparseUtils.jl.This package implements the type SparseMatrixCOO, several functions for type SparseMatrixCSC, and a few functions for AbstractSparseMatrix. Submodules NZVal and IJV contain the generic algorithms, which are wrapped by the sparse matrix types.Several functions are not yet documented."
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
    "location": "index.html#SparseUtils.SparseMatrixCOO",
    "page": "SparseUtils",
    "title": "SparseUtils.SparseMatrixCOO",
    "category": "type",
    "text": "SparseMatrixCOO{Tv,Ti<:Integer} <: AbstractSparseMatrix{Tv,Ti}\n\nMatrix type for storing sparse matrices in the Coordinate list form (COO) format.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseArrays.sparse",
    "page": "SparseUtils",
    "title": "SparseArrays.sparse",
    "category": "function",
    "text": "sparse(I, J, V, m, n; sparsetype=SparseMatrixCOO)\n\nCreate a sparse matrix of type SparseMatrixCOO\n\n\n\n\n\n"
},

{
    "location": "index.html#Base.issorted",
    "page": "SparseUtils",
    "title": "Base.issorted",
    "category": "function",
    "text": "issorted(S::SparseMatrixCOO)\n\nReturn true if the data structures in S storing non-zero indices are canonically sorted. If the data structures are not sorted, then functions taking an SparseMatrixCOO argument will, in general, give incorrect results.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseArrays.spzeros",
    "page": "SparseUtils",
    "title": "SparseArrays.spzeros",
    "category": "function",
    "text": "spzeros([type,] m,n; sparsetype=SparseMatrixCOO)\n\nCreate a sparse matrix of size mxn. If sparsetype is ommited, a matrix of type SparseMatrixCSC will be created.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.prunecols",
    "page": "SparseUtils",
    "title": "SparseUtils.prunecols",
    "category": "function",
    "text": "prunecols(S::SparseMatrixCOO, min_entries; renumber=true)\n\nRemove all columns from S that have fewer than min_entries non-zero entries. If renumber is true, renumber the remaining columns consecutively starting from 1.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.prunerows",
    "page": "SparseUtils",
    "title": "SparseUtils.prunerows",
    "category": "function",
    "text": "prunerows(S::SparseMatrixCOO, min_entries; renumber=true)\n\nRemove all rows from S that have fewer than min_entries non-zero entries. If renumber is true, renumber the remaining rows consecutively starting from 1.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseMatrixCOO-1",
    "page": "SparseUtils",
    "title": "SparseMatrixCOO",
    "category": "section",
    "text": "SparseUtils.SparseMatrixCOO\nSparseUtils.sparse\nSparseUtils.issorted\nSparseUtils.spzeros\nSparseUtils.prunecols\nSparseUtils.prunerows"
},

{
    "location": "index.html#SparseArrays.nnz-Tuple{SparseArrays.SparseMatrixCSC,Integer}",
    "page": "SparseUtils",
    "title": "SparseArrays.nnz",
    "category": "method",
    "text": "nnz(sp::SparseMatrixCSC, colnum::Integer)\n\nReturn the number of structurally non-zero elements in the colnumth column of sp.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.nnzcounts",
    "page": "SparseUtils",
    "title": "SparseUtils.nnzcounts",
    "category": "function",
    "text": "nnzcounts(sp::SparseMatrixCSC; rev=true, byvalue=true)\n\nReturn a sorted countmap of the number of structural non-zeros in the columns of sp. The keys controlling the sorting are passed to StatsBase.sort.\n\nExamples\n\nIn this example, there are 132781 columns with one non-zero element, 32096 columns with two, etc.\n\njulia> conns = nnzcounts(sp);\n\njulia> collect(Iterators.take(conns,3))\n3-element Array{Pair{Int64,Int64},1}:\n 1 => 132781\n 2 => 32096\n 3 => 17478\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.c_to_julia_index!",
    "page": "SparseUtils",
    "title": "SparseUtils.c_to_julia_index!",
    "category": "function",
    "text": "c_to_julia_index!(colptr, rowval, nzval)::SparseArrays.SparseMatrixCSC\n\nConvert a sparse matrix with zero-based indices to a SparseMatrixCSC. colptr and rowval are altered in place. nzval is neither copied nor altered.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.c_to_julia_index",
    "page": "SparseUtils",
    "title": "SparseUtils.c_to_julia_index",
    "category": "function",
    "text": "c_to_julia_index(colptr, rowval, nzval)::SparseArrays.SparseMatrixCSC\n\nConvert a sparse matrix with zero-based indices to a SparseMatrixCSC. colptr and rowval are not altered, but rather copied. nzval is neither copied, nor altered.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseMatrixCSC-1",
    "page": "SparseUtils",
    "title": "SparseMatrixCSC",
    "category": "section",
    "text": "SparseUtils.nnz(sp::SparseMatrixCSC, colnum::Integer)\nSparseUtils.nnzcounts\nSparseUtils.c_to_julia_index!\nSparseUtils.c_to_julia_index<!– SparseUtils.transpose –>"
},

{
    "location": "index.html#SparseUtils.density",
    "page": "SparseUtils",
    "title": "SparseUtils.density",
    "category": "function",
    "text": "density(sp::AbstractSparseMatrix)\n\nCompute the density of sp, the fraction of elements that are non-zero This assumes that there are no stored values equal to zero.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.summarystats",
    "page": "SparseUtils",
    "title": "SparseUtils.summarystats",
    "category": "function",
    "text": "summarystats(sp::SparseMatrixCSC)\n\nPrint some statistics for sp. The line \"number of non-zeros\" counts structural non-zeros.\n\n\n\n\n\n"
},

{
    "location": "index.html#AbstractSparseMatrix-1",
    "page": "SparseUtils",
    "title": "AbstractSparseMatrix",
    "category": "section",
    "text": "SparseUtils.density\nSparseUtils.summarystats"
},

{
    "location": "index.html#SparseUtils.NZVal",
    "page": "SparseUtils",
    "title": "SparseUtils.NZVal",
    "category": "module",
    "text": "module SparseUtils.NZVal\n\nFor sparse matrices, algorithms that depend only on nzval, m, and n.\n\nThis code is generic, for types in Core.\n\nThe symbols in NZVal should not be exported, not even for convenience. They conflict with functions in Base.\n\n\n\n\n\n"
},

{
    "location": "index.html#NZVal-1",
    "page": "SparseUtils",
    "title": "NZVal",
    "category": "section",
    "text": "SparseUtils.NZVal"
},

{
    "location": "index.html#SparseUtils.IJV",
    "page": "SparseUtils",
    "title": "SparseUtils.IJV",
    "category": "module",
    "text": "module IJV\n\nAlgorithms for sparse matrices stored in COO form.\n\nIJV does not depend on any types not in Core. Do not import any symbol from IJV into another module. The methods in IJV are incompatible with functions in Base with the same (unqualified) name. Instead, use IJV like this\n\n#Example\n\nBase.getindex(S::MySparseMatrixCOO, i::Integer, j::Integer) = IJV.getindex(S.I, S.J, S.V, i, j)\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.IJV.issorted",
    "page": "SparseUtils",
    "title": "SparseUtils.IJV.issorted",
    "category": "function",
    "text": "issorted(I::AbstractVector, J::AbstractVector)\n\nReturn true if the COO indices I and J are sorted canonically.\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.IJV.findindex",
    "page": "SparseUtils",
    "title": "SparseUtils.IJV.findindex",
    "category": "function",
    "text": "findindex(I::AbstractArray, J::AbstractArray, i::Integer, j::Integer)\n\nReturn the common linear index into I and J (and V) representing coordinates (i, j), or nothing if there is no value stored for (i, j).\n\n\n\n\n\n"
},

{
    "location": "index.html#SparseUtils.IJV.getindex",
    "page": "SparseUtils",
    "title": "SparseUtils.IJV.getindex",
    "category": "function",
    "text": "getindex(I, J, V, i::Integer, j::Integer; zeroel=zero(eltype(V)))\n\nReturn the element at coordinates (i, j) in the matrix represented by  I, J, V, if it exists, and zeroel otherwise.\n\n\n\n\n\n"
},

{
    "location": "index.html#IJV-1",
    "page": "SparseUtils",
    "title": "IJV",
    "category": "section",
    "text": "SparseUtils.IJV\nSparseUtils.IJV.issorted\nSparseUtils.IJV.findindex\nSparseUtils.IJV.getindex"
},

]}
