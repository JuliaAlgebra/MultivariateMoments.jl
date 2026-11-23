using Test, LinearAlgebra, MultivariateMoments

@testset "low_rank_ldlt" begin
    v = [1, -1]
    M = v * v'
    ldlt = low_rank_ldlt(M, SVDLDLT(), 1e-10)
    @test ldlt.L ≈ -normalize(v)
    ldlt = low_rank_ldlt(M, ShiftCholeskyLDLT(1e-10), 1e-9)
    @test ldlt.L ≈ v
    v = [1, -1 + im]
    M = v * v'
    ldlt = low_rank_ldlt(M, SVDLDLT(), 1e-10)
    @test ldlt.L ≈ -normalize(v)
end

@testset "Low-rank accuracy" begin
    check = LeadingRelativeRankTol(1e-8)
    for T in (Float32, Float64, BigFloat)
        L = ones(T, 2, 1)
        @test accuracy(LowRankLDLT(L, T[1, 0], check)) == 2eps(T)
        @test accuracy(LowRankLDLT(L, T[1, 1e-4], check)) ≈ T(1e-4)
    end
    # Gap right after the rank
    @test accuracy([1.0, 1e-3, 1e-9], 2, check) ≈ 1e-6
    # A later gap gives the same tolerance as a spectrum without that gap.
    @test accuracy([1.0, 1e-3, 1e-4, 1e-9], 2, check) ≈ 1e-4
    # No gap
    @test accuracy([1.0, 1e-3, 5e-4, 2e-4], 2, check) ≈ 5e-4
    # The noise `1e-20` is smaller than the error `3eps` of the singular values
    @test accuracy([1.0, 0.5, 1e-20], 2, check) ≈ 3eps() / 0.5
end

@testset "Configurable spectral gap" begin
    σ = [1.0, 0.125, 0.015625, 1e-6, 1e-10]
    check = FixedRank(2)
    @test accuracy(σ, 2, GapBasedAccuracy(check)) == accuracy(σ, 2, check)
    # The comparison at the threshold is strict (σ[3] / σ[2] = 0.125).
    @test accuracy(σ, 2, GapBasedAccuracy(check; gap_tol = 0.125)) == σ[3]
    @test accuracy(σ, 2, GapBasedAccuracy(check; gap_tol = 0.25)) == σ[3] / σ[2]
    later = GapBasedAccuracy(check; search_later = true)
    # Take the first qualifying gap, even when another one follows.
    @test accuracy(σ, 2, later) == σ[4] / σ[2]
    @test accuracy(
        σ,
        2,
        GapBasedAccuracy(check; gap_tol = 0, search_later = true),
    ) == σ[3]
    @test accuracy([1.0, 0.125, 0.015625, 0.0], 2, later) == 4eps() / 0.125
    @test accuracy(σ, length(σ), later) == accuracy(σ, length(σ), check)
    @test iszero(accuracy(σ, 0, later))
    chol = low_rank_ldlt(Diagonal(σ), SVDLDLT(), later)
    @test size(chol.L, 2) == 2
    @test accuracy(chol) == σ[4] / σ[2]
    @test doubt(σ, GapBasedAccuracy(LargestDifferentialRank())) ==
          doubt(σ, LargestDifferentialRank())
end

struct HardcodedRanks <: RankCheck
    r::Vector{Int}
end

function MultivariateMoments.rank_from_singular_values(σ, r::HardcodedRanks)
    return r.r[length(σ)]
end

@testset "Decreasing rank" begin
    r = RankDependence(zeros(4, 4), HardcodedRanks([1, 0, 1, 2]))
    @test !is_dependent!(r, 1)
    @test is_dependent!(r, 2)
    @test is_dependent!(r, 3)
    @test is_dependent!(r, 4)
end
