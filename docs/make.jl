using Documenter, SparseUtils

makedocs(
    debug = true,
    strict = false,
    doctest = false,
    format = :html,
    sitename = "SparseUtils.jl",
    modules = [SparseUtils],
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo = "github.com/jlapeyre/SparseUtils.jl.git",
    target = "build",
    julia  = "0.7",
    deps = nothing,
    make = nothing
)
