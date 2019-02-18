using Documenter, MultivariateMoments

makedocs(
    sitename = "MultivariateMoments",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
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
)
