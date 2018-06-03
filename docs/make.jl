using Documenter, MultivariateMoments

makedocs(
    format = :html,
    sitename = "MultivariateMoments",
    pages = [
        "Index" => "index.md"
        "Moments and expectation" => "moments.md"
        "Atoms extraction" => "atoms.md"
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [MultivariateMoments]
)

deploydocs(
    repo   = "github.com/JuliaAlgebra/MultivariateMoments.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
