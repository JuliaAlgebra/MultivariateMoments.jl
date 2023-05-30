using SemialgebraicSets

abstract type AbstractMomentMatrix{T,B<:MB.AbstractPolynomialBasis} <:
              AbstractMeasureLike{T} end

"""
    mutable struct MomentMatrix{T, B<:MultivariateBases.AbstractPolynomialBasis, MT<:AbstractMatrix{T}} <: AbstractMeasureLike{T}
        Q::MT
        basis::B
        support::Union{Nothing, AlgebraicSet}
    end

Measure ``\\nu`` represented by the moments of the monomial matrix ``x x^\\top`` in the symmetric matrix `Q`.
The set of points that are zeros of all the polynomials ``p`` such that ``\\mathbb{E}_{\\nu}[p] = 0`` is stored in the field `support` when it is computed.
"""
mutable struct MomentMatrix{
    T,
    B<:MB.AbstractPolynomialBasis,
    MT<:AbstractMatrix{T},
} <: AbstractMomentMatrix{T,B}
    Q::MT
    basis::B
    support::Union{Nothing,AbstractAlgebraicSet}
end
function MomentMatrix{T,B,MT}(
    Q::MT,
    basis::MB.AbstractPolynomialBasis,
) where {T,B,MT}
    return MomentMatrix{T,B,MT}(Q, basis, nothing)
end
function MomentMatrix{T,B}(
    Q::AbstractMatrix{T},
    basis::MB.AbstractPolynomialBasis,
) where {T,B}
    return MomentMatrix{T,B,typeof(Q)}(Q, basis)
end
function MomentMatrix(
    Q::SymMatrix{T},
    basis::MB.AbstractPolynomialBasis,
) where {T}
    return MomentMatrix{T,typeof(basis)}(Q, basis)
end

MP.variables(μ::MomentMatrix) = MP.variables(μ.basis)
MP.nvariables(μ::MomentMatrix) = MP.nvariables(μ.basis)

function MomentMatrix{T}(
    f::Function,
    basis::MB.AbstractPolynomialBasis,
    σ = 1:length(basis),
) where {T}
    return MomentMatrix(
        vectorized_symmetric_matrix(T, f, length(basis), σ),
        basis,
    )
end
function MomentMatrix{T}(f::Function, monos::AbstractVector) where {T}
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return MomentMatrix{T}(f, MB.MonomialBasis(sorted_monos), σ)
end

"""
    moment_matrix(μ::Measure, x)

Creates a matrix the moment matrix for the moment matrix  ``x x^\\top`` using the moments of `μ`.
"""
function moment_matrix(μ::Measure{T}, X) where {T}
    return MomentMatrix{T}((i, j) -> moment_value(μ, X[i] * X[j]), X)
end

function MomentMatrix(
    Q::AbstractMatrix{T},
    basis::MB.AbstractPolynomialBasis,
    σ,
) where {T}
    return MomentMatrix{T}((i, j) -> Q[σ[i], σ[j]], basis)
end
function MomentMatrix(Q::AbstractMatrix, monos::AbstractVector)
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return MomentMatrix(Q, MB.MonomialBasis(sorted_monos), σ)
end
moment_matrix(Q::AbstractMatrix, monos) = MomentMatrix(Q, monos)

value_matrix(μ::MomentMatrix) = Matrix(μ.Q)

function vectorized_basis(ν::MomentMatrix{T,<:MB.MonomialBasis}) where {T}
    monos = ν.basis.monomials
    n = length(monos)
    return MB.MonomialBasis([monos[i] * monos[j] for i in 1:n for j in 1:i])
end

function measure(ν::MomentMatrix; kws...)
    n = length(ν.basis)
    return measure(ν.Q.Q, vectorized_basis(ν); kws...)
end

struct SparseMomentMatrix{T,B<:MB.AbstractPolynomialBasis,MT} <:
       AbstractMomentMatrix{T,B}
    sub_moment_matrices::Vector{MomentMatrix{T,B,MT}}
end
