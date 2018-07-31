#__precompile__()
module SparseUtils

import SparseArrays
import SparseArrays: SparseMatrixCSC, nnz
import LinearAlgebra
import DataUtils
import StatsBase
import Printf

const Callable = Union{Function, Type}

export
    c_to_julia_index,
    c_to_julia_index!,
    density,
    nnzcounts,
    summarystats,
    nrows,
    ncols,
    prunecols!,
    prunecols,
    prunerows,
    renumbercols,
    renumberrows,
    renumberrowscols,
    hasemptycols,
    hasemptyrows

include("findrepeated.jl")
include("nzval.jl")
include("ijv.jl")
include("csc.jl")
include("coo.jl")
include("coorenumber.jl")
include("abstract.jl")

import .NZVal, .IJV

end # module SparseUtils

#  LocalWords:  nnzcounts countmap sp
