using Documenter, Literate, MultivariateMoments
import DocumenterCitations
# Needed because we include docstrings for `monomial`, `variables` and `monomials`
import MultivariatePolynomials

DocMeta.setdocmeta!(
    MultivariateMoments,
    :DocTestSetup,
    :(using MultivariateMoments);
    recursive = true,
)

const _TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorials")
const _OUTPUT_DIR = joinpath(@__DIR__, "src", "generated")

for filename in readdir(_TUTORIAL_DIR)
    path = joinpath(_TUTORIAL_DIR, filename)
    Literate.markdown(path, _OUTPUT_DIR; documenter = true)
    Literate.notebook(path, _OUTPUT_DIR; documenter = true)
    Literate.script(path, _OUTPUT_DIR; documenter = true)
end

const _TUTORIALS = [
    joinpath("generated", file) for
    file in sort(readdir(_OUTPUT_DIR)) if endswith(file, ".md")
]

makedocs(
    sitename = "MultivariateMoments",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/tables.css"],
    ),
    pages = [
        "Index" => "index.md"
        "Moments and expectation" => "moments.md"
        "Atoms extraction" => "atoms.md"
        "Tutorials" => _TUTORIALS
        "Bibliography" => "bibliography.md"
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base and MultivariatePolynomials that
    # we overwrite.
    modules = [MultivariateMoments],
    plugins = [
        DocumenterCitations.CitationBibliography(
            joinpath(@__DIR__, "src", "references.bib");
            style = :authoryear,
        ),
    ],
)

deploydocs(
    repo = "github.com/JuliaAlgebra/MultivariateMoments.jl.git",
    push_preview = true,
)
