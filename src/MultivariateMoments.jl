module MultivariateMoments

using LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials

import MultivariateBases
const MB = MultivariateBases

export AbstractMeasureLike, AbstractMomentLike
abstract type AbstractMeasureLike{T} end
abstract type AbstractMomentLike{T} <: AbstractMeasureLike{T} end

export AbstractMeasure, AbstractMoment
abstract type AbstractMoment{T} <: AbstractMomentLike{T} end
abstract type AbstractMeasure{T} <: AbstractMeasureLike{T} end

include("moment.jl")
include("measure.jl")
include("expectation.jl")
include("symmatrix.jl")
include("moment_matrix.jl")
include("atomic.jl")
include("extract.jl")

end # module
