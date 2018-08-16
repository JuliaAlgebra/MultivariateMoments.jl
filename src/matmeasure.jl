export SymMatrix, MatMeasure, getmat, matmeasure

using SemialgebraicSets

struct SymMatrix{T} <: AbstractMatrix{T}
    Q::Vector{T}
    n
end

# i <= j
#function trimap(i, j, n) # It was used when Q was the lower triangular part
#    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
#end

# j <= i
function trimap(i, j)
    div((i-1)*i, 2) + j
end

function trimat(::Type{T}, f, n, σ) where {T}
    Q = Vector{T}(undef, trimap(n, n))
    for i in 1:n
        for j in 1:i
            Q[trimap(i, j)] = f(σ[i], σ[j])
        end
    end
    SymMatrix{T}(Q, n)
end

Base.size(Q::SymMatrix) = (Q.n, Q.n)

function Base.getindex(Q::SymMatrix, i::Integer, j::Integer)
    Q.Q[trimap(max(i, j), min(i, j))]
end
Base.getindex(Q::SymMatrix, I::Tuple) = Q[I...]
Base.getindex(Q::SymMatrix, I::CartesianIndex) = Q[I.I]

"""
    mutable struct MatMeasure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractMeasureLike{T}
        Q::SymMatrix{T}
        x::MVT
        support::Union{Nothing, AlgebraicSet}
    end

Measure ``\\nu`` represented by the moments of the monomial matrix ``x x^\\top`` in the symmetric matrix `Q`.
The set of points that are zeros of all the polynomials ``p`` such that ``\\mathbb{E}_{\\nu}[p] = 0`` is stored in the field `support` when it is computed.
"""
mutable struct MatMeasure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractMeasureLike{T}
    Q::SymMatrix{T}
    x::MVT
    support::Union{Nothing, AlgebraicSet}
end
MatMeasure{T, MT, MVT}(Q::SymMatrix{T}, x::MVT) where {T, MT, MVT} = MatMeasure{T, MT, MVT}(Q, x, nothing)

MP.variables(μ::MatMeasure) = variables(μ.x)
MP.nvariables(μ::MatMeasure) = nvariables(μ.x)

function MatMeasure{T}(f::Function, x::AbstractVector{MT}, σ) where {T, MT<:AbstractMonomial}
    MatMeasure{T, MT, monovectype(x)}(trimat(T, f, length(x), σ), x)
end
MatMeasure{T}(f::Function, x::AbstractVector, σ) where T = MatMeasure{T}(f, monomial.(x), σ)
function MatMeasure{T}(f::Function, x::AbstractVector) where T
    σ, X = sortmonovec(x)
    MatMeasure{T}(f, X, σ)
end
MatMeasure(f::Function, x) = MatMeasure{Base.promote_op(f, Int, Int)}(f, x)
matmeasure(f::Function, x) = MatMeasure(f, x)

"""
    matmeasure(μ::Measure, x)

Creates a matrix the moment matrix for the moment matrix  ``x x^\\top`` using the moments of `μ`.
"""
function matmeasure(μ::Measure{T}, X) where T
    function getmom(i, j)
        x = X[i] * X[j]
        for m in moments(μ)
            if monomial(m) == x
                return value(m)
            end
        end
        throw(ArgumentError("μ does not have the moment $(x)"))
    end
    MatMeasure{T}(getmom, X)
end

function MatMeasure(Q::AbstractMatrix{T}, x, σ) where T
    MatMeasure{T}((i,j) -> Q[σ[i], σ[j]], x)
end
function MatMeasure(Q::AbstractMatrix, x)
    σ, X = sortmonovec(x)
    MatMeasure(Q, X, σ)
end
matmeasure(Q::AbstractMatrix, x) = MatMeasure(Q, x)

function getmat(μ::MatMeasure)
    Matrix(μ.Q)
end

function measure(ν::MatMeasure)
    n = length(ν.x)
    measure(ν.Q.Q, [ν.x[i] * ν.x[j] for i in 1:n for j in 1:i])
end
