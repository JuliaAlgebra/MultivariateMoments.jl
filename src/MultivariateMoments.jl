__precompile__()

module MultivariateMoments

using MultivariatePolynomials
const MP = MultivariatePolynomials

include("measure.jl")
include("exp.jl")
include("matmeasure.jl")

end # module
