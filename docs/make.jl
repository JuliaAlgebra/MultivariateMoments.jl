using Documenter, MultivariateMoments
# Needed because we include docstrings for `monomial`, `variables` and `monomials`
import MultivariatePolynomials

DocMeta.setdocmeta!(MultivariateMoments, :DocTestSetup, :(using MultivariateMoments); recursive=true)

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
    # this module for functions define in Base and MultivariatePolynomials that
    # we overwrite.
    modules = [MultivariateMoments],
)

deploydocs(
    repo   = "github.com/JuliaAlgebra/MultivariateMoments.jl.git",
)
