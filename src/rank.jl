export LowRankChol, ShiftCholesky, SVDCholesky
export FixedRank, AbsoluteRankTol, LeadingRelativeRankTol, DifferentialRankTol, LargestDifferentialRank

abstract type RankCheck end

struct FixedRank <: RankCheck
    r::Int
end
rank_from_singular_values(σ, check::FixedRank) = check.r

function _findfirst(f, σ)
    return something(findfirst(f, eachindex(σ)), length(σ) + 1) - 1
end

struct AbsoluteRankTol{T} <: RankCheck
    tol::T
end
rank_from_singular_values(σ, check::AbsoluteRankTol) = _findfirst(i -> σ[i] ≤ check.tol, σ)

struct LeadingRelativeRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::LeadingRelativeRankTol)
    return _findfirst(σ) do i
        return σ[i] ≤ check.tol * first(σ)
    end
end

struct DifferentialRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::DifferentialRankTol)
    return _findfirst(σ) do i
        if i == 1
            return false
        else
            return σ[i] ≤ check.tol * σ[i - 1]
        end
    end
end

struct LargestDifferentialRank <: RankCheck end
function rank_from_singular_values(σ, ::LargestDifferentialRank)
    return argmax(i -> σ[i] / σ[i + 1], 1:(length(σ) - 1))
end


"""
    LowRankChol

Method for computing a ``r \\times n`` matrix `U` of a ``n \\times n`` rank ``r`` psd matrix `Q` such that `Q = U'U`.
"""
abstract type LowRankChol end

"""
    ShiftCholesky <: LowRankChol

Shift the matrix by `shift` times the identity matrix before cholesky.
"""
struct ShiftCholesky{T} <: LowRankChol
    shift::T
end

"""
    SVDCholesky <: LowRankChol

Use SVD decomposition.
"""
struct SVDCholesky <: LowRankChol end

"""
    MultivariateMoments.low_rank_cholesky(Q::AbstractMatrix, dec::LowRankChol, ranktol)

Returns a ``r \\times n`` matrix ``U`` of a ``n \\times n`` rank ``r`` positive semidefinite matrix ``Q`` such that ``Q = U^\\top U``.
The rank of ``Q`` is the number of singular values larger than `ranktol```{} \\cdot \\sigma_1`` where ``\\sigma_1`` is the largest singular value.
"""
function low_rank_cholesky end

"""
    struct LowRankCholesky{T}
        U::Matrix{T}
        singular_values::Vector{T}
    end

Low-Rank cholesky decomposition `U` of size `(r, n)` of a `n x n` matrix
with singular values `singular_values[1] > ... > singular_values[n]`.
"""
struct LowRankCholesky{T}
    U::Matrix{T}
    singular_values::Vector{T}
end

"""
    spectral_norm(chol::LowRankCholesky)

Return the spectral norm of the matrix `chol.U' * chol.U`.
"""
spectral_norm(chol::LowRankCholesky) = first(chol.singular_values)

"""
    recommended_rtol(chol::LowRankCholesky)

Return the ratio `rtol = σ[r+1]/σ[1]` where `σ` is the vector of singular
values of the matrix for which `chol` is the rank-`r` Cholesky decomposition.
This is a good relative tolerance to use with the matrix as `σ[r+1]` is the
first singular value that was discarded.
"""
function recommended_rtol(chol::LowRankCholesky)
    r = size(chol.U, 1)
    first_dropped = min(length(chol.singular_values), r + 1)
    # TODO If `r = length(chol.singular_values)`, we would like to take `rank_check.tol`
    # but some don't have any rank check, we could still use it when there is a tol
    # Also, dividing by `nM` is maybe not the best if we don't use `LeadingRelativeRankTol`
    return chol.singular_values[first_dropped] / spectral_norm(chol)
end

function low_rank_cholesky(M::AbstractMatrix, dec::ShiftCholesky, rank_check::RankCheck)
    m = LinearAlgebra.checksquare(M)
    U = cholesky(M + dec.shift * I).U
    σs = map(i -> (U[i, i])^2, 1:m)
    J = sortperm(σs, rev=true)
    σ_sorted = σs[J]
    r = rank_from_singular_values(σ_sorted, rank_check)
    return LowRankCholesky(U[J[1:r], :], σ_sorted)
end

function low_rank_cholesky(M::AbstractMatrix, ::SVDCholesky, rank_check::RankCheck)
    F = svd(M)
    r = rank_from_singular_values(F.S, rank_check)
    return LowRankCholesky(Diagonal(sqrt.(F.S[1:r])) * F.U[:, 1:r]', F.S)
end
