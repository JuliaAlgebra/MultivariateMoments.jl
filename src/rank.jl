abstract type RankCheck end

"""
    rank_from_singular_values(σ, check::RankCheck)

Return the rank of a matrix with singular values `σ` (in decreasing order)
using `check`.
"""
function rank_from_singular_values end

function accuracy(σ::Vector{T}, r, check::RankCheck) where {T}
    if iszero(r)
        return zero(T)
    elseif r == lastindex(σ)
        return default_accuracy(σ, check)
    else
        return σ[r+1] / σ[1]
    end
end

struct FixedRank <: RankCheck
    r::Int
end
rank_from_singular_values(_, check::FixedRank) = check.r
default_accuracy(::Vector{T}, ::FixedRank) where {T} = Base.rtoldefault(T)

function _findfirst(f, σ)
    return something(findfirst(f, eachindex(σ)), lastindex(σ) + 1) - 1
end

struct AbsoluteRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::AbsoluteRankTol)
    return _findfirst(i -> σ[i] ≤ check.tol, σ)
end
default_accuracy(σ::Vector, check::AbsoluteRankTol) = check.tol / last(σ)

struct LeadingRelativeRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::LeadingRelativeRankTol)
    return _findfirst(σ) do i
        return σ[i] ≤ check.tol * first(σ)
    end
end
function default_accuracy(σ::Vector, check::LeadingRelativeRankTol)
    return check.tol * first(σ) / last(σ)
end

struct DifferentialRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::DifferentialRankTol)
    return _findfirst(σ) do i
        if i == 1
            return false
        else
            return σ[i] ≤ check.tol * σ[i-1]
        end
    end
end
default_accuracy(::Vector, check::DifferentialRankTol) = check.tol

struct LargestDifferentialRank <: RankCheck end
function rank_from_singular_values(σ, ::LargestDifferentialRank)
    return argmax(i -> σ[i] / σ[i+1], 1:(length(σ)-1))
end
function default_accuracy(::Vector{T}, ::LargestDifferentialRank) where {T}
    return Base.rtoldefault(T)
end

"""
    LowRankLDLTAlgorithm

Method for computing a ``r \\times n`` matrix `U` of a ``n \\times n`` rank ``r`` psd matrix `Q` such that `Q = U'U`.
"""
abstract type LowRankLDLTAlgorithm end

"""
    ShiftCholeskyLDLT <: LowRankLDLTAlgorithm

Shift the matrix by `shift` times the identity matrix before cholesky.
"""
struct ShiftCholeskyLDLT{T} <: LowRankLDLTAlgorithm
    shift::T
end

"""
    SVDLDLT <: LowRankLDLTAlgorithm

Use SVD decomposition.
"""
struct SVDLDLT <: LowRankLDLTAlgorithm end

"""
    MultivariateMoments.low_rank_ldlt(Q::AbstractMatrix, dec::LowRankLDLTAlgorithm, ranktol)

Returns a ``r \\times n`` matrix ``U`` of a ``n \\times n`` rank ``r`` positive semidefinite matrix ``Q`` such that ``Q = U^\\top U``.
The rank of ``Q`` is the number of singular values larger than `ranktol```{} \\cdot \\sigma_1`` where ``\\sigma_1`` is the largest singular value.
"""
function low_rank_ldlt end

"""
    struct LowRankLDLT{T}
        L::Matrix{T}
        singular_values::Vector{T}
    end

Low-Rank cholesky decomposition `L * Diagonal(singular_values) * L'` of size
`(n, r)` of a `n x n` matrix with singular values
`singular_values[1] > ... > singular_values[n]`.
The rank was chosen given `singular_values` using `rank_check` via
the [`rank_from_singular_values`](@ref) function.
"""
struct LowRankLDLT{T,C<:RankCheck}
    L::Matrix{T}
    singular_values::Vector{T}
    rank_check::C
end

"""
    accuracy(chol::LowRankLDLT)

Return the ratio `rtol = σ[r+1]/σ[1]` where `σ` is the vector of singular
values of the matrix for which `chol` is the rank-`r` Cholesky decomposition.
This is a good relative tolerance to use with the matrix as `σ[r+1]` is the
first singular value that was discarded.
"""
function accuracy(chol::LowRankLDLT)
    return accuracy(chol.singular_values, size(chol.L, 2), chol.rank_check)
end

function low_rank_ldlt(
    M::AbstractMatrix,
    dec::ShiftCholeskyLDLT,
    rank_check::RankCheck,
)
    m = LinearAlgebra.checksquare(M)
    U = cholesky(M + dec.shift * I).U
    σs = map(i -> (U[i, i])^2, 1:m)
    J = sortperm(σs, rev = true)
    σ_sorted = σs[J]
    r = rank_from_singular_values(σ_sorted, rank_check)
    return LowRankLDLT(
        U[J[1:r], :]' * Diagonal(inv.(sqrt.(σ_sorted[1:r]))),
        σ_sorted,
        rank_check,
    )
end

function low_rank_ldlt(M::AbstractMatrix, ::SVDLDLT, rank_check::RankCheck)
    F = svd(M)
    r = rank_from_singular_values(F.S, rank_check)
    return LowRankLDLT(F.U[:, 1:r], F.S, rank_check)
end
