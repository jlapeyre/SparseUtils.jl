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
    julia  = "nightly",
    deps = nothing,
    make = nothing
)
