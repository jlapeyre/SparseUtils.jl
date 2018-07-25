__precompile__()

module SparseUtils

import SparseArrays
import SparseArrays: SparseMatrixCSC, nnz
import DataUtils
import StatsBase
import Printf

export
    c_to_julia_index,
    c_to_julia_index!,
    sparsity,
    nnzcounts,
    sparse_stats,
    nrows,
    ncols,
    prunecols!,
    dropzerocols,
    dropzerorows,
    dropzerorowscols,
    hasemptycols,
    hasemptyrows

include("csc.jl")

end # module SparseUtils

#  LocalWords:  nnzcounts countmap sp
