export SymMatrix, MomentMatrix, getmat, moment_matrix, symmetric_setindex!

using SemialgebraicSets

"""
    mutable struct MomentMatrix{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractMeasureLike{T}
        Q::SymMatrix{T}
        x::MVT
        support::Union{Nothing, AlgebraicSet}
    end

Measure ``\\nu`` represented by the moments of the monomial matrix ``x x^\\top`` in the symmetric matrix `Q`.
The set of points that are zeros of all the polynomials ``p`` such that ``\\mathbb{E}_{\\nu}[p] = 0`` is stored in the field `support` when it is computed.
"""
mutable struct MomentMatrix{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractMeasureLike{T}
    Q::SymMatrix{T}
    x::MVT
    support::Union{Nothing, AlgebraicSet}
end
MomentMatrix{T, MT, MVT}(Q::SymMatrix{T}, x::MVT) where {T, MT, MVT} = MomentMatrix{T, MT, MVT}(Q, x, nothing)

MP.variables(μ::MomentMatrix) = variables(μ.x)
MP.nvariables(μ::MomentMatrix) = nvariables(μ.x)

function MomentMatrix{T}(f::Function, x::AbstractVector{MT}, σ) where {T, MT<:AbstractMonomial}
    MomentMatrix{T, MT, monovectype(x)}(trimat(T, f, length(x), σ), x)
end
MomentMatrix{T}(f::Function, x::AbstractVector, σ) where T = MomentMatrix{T}(f, monomial.(x), σ)
function MomentMatrix{T}(f::Function, x::AbstractVector) where T
    σ, X = sortmonovec(x)
    MomentMatrix{T}(f, X, σ)
end
MomentMatrix(f::Function, x) = MomentMatrix{Base.promote_op(f, Int, Int)}(f, x)
moment_matrix(f::Function, x) = MomentMatrix(f, x)

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
    MomentMatrix{T}(getmom, X)
end

function MomentMatrix(Q::AbstractMatrix{T}, x, σ) where T
    MomentMatrix{T}((i,j) -> Q[σ[i], σ[j]], x)
end
function MomentMatrix(Q::AbstractMatrix, x)
    σ, X = sortmonovec(x)
    MomentMatrix(Q, X, σ)
end
moment_matrix(Q::AbstractMatrix, x) = MomentMatrix(Q, x)

function getmat(μ::MomentMatrix)
    Matrix(μ.Q)
end

function measure(ν::MomentMatrix)
    n = length(ν.x)
    measure(ν.Q.Q, [ν.x[i] * ν.x[j] for i in 1:n for j in 1:i])
end
