module MultivariateMoments

using LinearAlgebra

import MutableArithmetics as MA

import MultivariatePolynomials as MP

import MultivariateBases as MB

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
include("hermitian_matrix.jl")
include("moment_matrix.jl")
include("atomic.jl")

include("rank.jl")
include("echelon.jl")
include("extract.jl")
include("shift.jl")
include("flat.jl")

end # module
