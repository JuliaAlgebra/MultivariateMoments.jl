__precompile__()

module MultivariateMoments

using MultivariatePolynomials
const MP = MultivariatePolynomials

export AbstractMeasureLike, AbstractMomentLike
abstract type AbstractMeasureLike{T} end
abstract type AbstractMomentLike{T} <: AbstractMeasureLike{T} end

abstract type AbstractMoment{T} <: AbstractMomentLike{T} end
abstract type AbstractMeasure{T} <: AbstractMeasureLike{T} end

import Base: *, +

include("moment.jl")
include("measure.jl")
include("exp.jl")
include("matmeasure.jl")

end # module
