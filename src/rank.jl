export LowRankChol, ShiftChol, SVDChol
export FixedRank,
    AbsoluteRankTol, LeadingRelativeRankTol, DifferentialRankTol, LargestDifferentialRank

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
            return σ[i] ≤ check.tol * σ[i-1]
        end
    end
end

struct LargestDifferentialRank <: RankCheck end
function rank_from_singular_values(σ, ::LargestDifferentialRank)
    return argmax(i -> σ[i] / σ[i+1], 1:(length(σ)-1))
end


"""
    LowRankChol

Method for computing a ``r \\times n`` matrix `U` of a ``n \\times n`` rank ``r`` psd matrix `Q` such that `Q = U'U`.
"""
abstract type LowRankChol end

"""
    ShiftChol <: LowRankChol

Shift the matrix by `shift` times the identity matrix before cholesky.
"""
struct ShiftChol{T} <: LowRankChol
    shift::T
end

"""
    SVDChol <: LowRankChol

Use SVD decomposition.
"""
struct SVDChol <: LowRankChol end

"""
    MultivariateMoments.lowrankchol(Q::AbstractMatrix, dec::LowRankChol, ranktol)

Returns a ``r \\times n`` matrix ``U`` of a ``n \\times n`` rank ``r`` positive semidefinite matrix ``Q`` such that ``Q = U^\\top U``.
The rank of ``Q`` is the number of singular values larger than `ranktol```{} \\cdot \\sigma_1`` where ``\\sigma_1`` is the largest singular value.
"""
function lowrankchol end

# TODO If `r = length(σ)`, we would like to take `rank_check.tol`
# but some don't have any rank check, we could still use it when there is a tol
# Also, dividing by `nM` is maybe not the best if we don't use `LeadingRelativeRankTol`
function _M(σ, r)
    nM = first(σ)
    cM = σ[min(length(σ), r + 1)] / nM
    return nM, cM
end

function lowrankchol(M::AbstractMatrix, dec::ShiftChol, rank_check::RankCheck)
    m = LinearAlgebra.checksquare(M)
    U = cholesky(M + dec.shift * I).U
    σs = map(i -> (U[i, i])^2, 1:m)
    J = sortperm(σs, rev = true)
    σ_sorted = σs[J]
    r = rank_from_singular_values(σ_sorted, rank_check)
    nM, cM = _M(σ_sorted, r)
    return nM, cM, U[J[1:r], :]
end

function lowrankchol(M::AbstractMatrix, ::SVDChol, rank_check::RankCheck)
    F = svd(M)
    r = rank_from_singular_values(F.S, rank_check)
    nM, cM = _M(F.S, r)
    return nM, cM, (F.U[:, 1:r] * Diagonal(sqrt.(F.S[1:r])))'
end
