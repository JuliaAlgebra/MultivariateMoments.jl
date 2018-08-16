__precompile__()

module MultivariateMoments

using Compat
using Compat.LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials

export AbstractMeasureLike, AbstractMomentLike
abstract type AbstractMeasureLike{T} end
abstract type AbstractMomentLike{T} <: AbstractMeasureLike{T} end

abstract type AbstractMoment{T} <: AbstractMomentLike{T} end
abstract type AbstractMeasure{T} <: AbstractMeasureLike{T} end

include("moment.jl")
include("measure.jl")
include("expectation.jl")
include("matmeasure.jl")
include("atomic.jl")
include("extract.jl")

end # module
