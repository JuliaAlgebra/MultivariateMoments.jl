export SymMatrix, AbstractMomentMatrix, MomentMatrix, SparseMomentMatrix
export getmat, moment_matrix, symmetric_setindex!

using SemialgebraicSets

abstract type AbstractMomentMatrix{T, B<:MB.AbstractPolynomialBasis} <: AbstractMeasureLike{T} end

"""
    mutable struct MomentMatrix{T, B<:MultivariateBases.AbstractPolynomialBasis} <: AbstractMeasureLike{T}
        Q::SymMatrix{T}
        basis::B
        support::Union{Nothing, AlgebraicSet}
    end

Measure ``\\nu`` represented by the moments of the monomial matrix ``x x^\\top`` in the symmetric matrix `Q`.
The set of points that are zeros of all the polynomials ``p`` such that ``\\mathbb{E}_{\\nu}[p] = 0`` is stored in the field `support` when it is computed.
"""
mutable struct MomentMatrix{T, B<:MB.AbstractPolynomialBasis} <: AbstractMomentMatrix{T, B}
    Q::SymMatrix{T}
    basis::B
    support::Union{Nothing, AlgebraicSet}
end
MomentMatrix{T, B}(Q::SymMatrix{T}, basis::MB.AbstractPolynomialBasis) where {T, B} = MomentMatrix{T, B}(Q, basis, nothing)
function MomentMatrix(Q::SymMatrix{T}, basis::MB.AbstractPolynomialBasis) where T
    return MomentMatrix{T, typeof(basis)}(Q, basis)
end

MP.variables(μ::MomentMatrix) = variables(μ.basis)
MP.nvariables(μ::MomentMatrix) = nvariables(μ.basis)

function MomentMatrix{T}(f::Function, basis::MB.AbstractPolynomialBasis, σ=1:length(basis)) where T
    return MomentMatrix(trimat(T, f, length(basis), σ), basis)
end
function MomentMatrix{T}(f::Function, monos::AbstractVector) where T
    σ, sorted_monos = sortmonovec(monos)
    return MomentMatrix{T}(f, MB.MonomialBasis(sorted_monos), σ)
end

"""
    moment_matrix(μ::Measure, x)

Creates a matrix the moment matrix for the moment matrix  ``x x^\\top`` using the moments of `μ`.
"""
function moment_matrix(μ::Measure{T}, X) where T
    function getmom(i, j)
        x = X[i] * X[j]
        for m in moments(μ)
            if monomial(m) == x
                return moment_value(m)
            end
        end
        throw(ArgumentError("μ does not have the moment $(x)"))
    end
    return MomentMatrix{T}(getmom, X)
end

function MomentMatrix(Q::AbstractMatrix{T}, basis::MB.AbstractPolynomialBasis, σ) where T
    return MomentMatrix{T}((i, j) -> Q[σ[i], σ[j]], basis)
end
function MomentMatrix(Q::AbstractMatrix, monos::AbstractVector)
    σ, sorted_monos = MP.sortmonovec(monos)
    return MomentMatrix(Q, MB.MonomialBasis(sorted_monos), σ)
end
moment_matrix(Q::AbstractMatrix, monos) = MomentMatrix(Q, monos)

function getmat(μ::MomentMatrix)
    Matrix(μ.Q)
end

function measure(ν::MomentMatrix{T, <:MB.MonomialBasis}) where T
    n = length(ν.basis)
    monos = ν.basis.elements
    measure(ν.Q.Q, [monos[i] * monos[j] for i in 1:n for j in 1:i])
end

struct SparseMomentMatrix{T, B <: MB.AbstractPolynomialBasis} <: AbstractMomentMatrix{T, B}
    sub_moment_matrices::Vector{MomentMatrix{T, B}}
end
