module MultivariateMomentsPositiveFactorizationsExt

using LinearAlgebra
import MultivariateMoments as MM
import PositiveFactorizations as PF

function MM.low_rank_ldlt(
    M::AbstractMatrix{<:Real},
    ::MM.AbsCholeskyLDLT,
    rank_check::MM.RankCheck,
)
    m = LinearAlgebra.checksquare(M)
    if iszero(m)
        T = float(eltype(M))
        return MM.LowRankLDLT(zeros(T, 0, 0), T[], rank_check)
    end
    F, signs = ldlt(PF.Positive, M)
    U = F.U
    # PositiveFactorizations replaces zero pivots by ones. They must not
    # contribute to the rank or to the reconstructed matrix.
    σs = [iszero(signs[i]) ? zero(eltype(U)) : U[i, i]^2 for i in 1:m]
    J = sortperm(σs; rev = true)
    σ_sorted = σs[J]
    r = MM.rank_from_singular_values(σ_sorted, rank_check)
    return MM.LowRankLDLT(
        U[J[1:r], :]' * Diagonal(inv.(sqrt.(σ_sorted[1:r]))),
        σ_sorted,
        rank_check,
    )
end

end
