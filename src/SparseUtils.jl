#__precompile__()
module SparseUtils

import SparseArrays
import SparseArrays: SparseMatrixCSC, nnz
import LinearAlgebra
import DataUtils
import StatsBase
import Printf

export
    c_to_julia_index,
    c_to_julia_index!,
    sparsity,
    nnzcounts,
    summarystats,
    nrows,
    ncols,
    prunecols!,
    renumbercols,
    renumberrows,
    renumberrowscols,
    hasemptycols,
    hasemptyrows

include("nzval.jl")
include("ijv.jl")
include("coo.jl")
include("csc.jl")
include("abstract.jl")


end # module SparseUtils

#  LocalWords:  nnzcounts countmap sp
