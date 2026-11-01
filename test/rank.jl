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
        @test accuracy(LowRankLDLT(L, T[1, 0], check)) == eps(T)
        @test accuracy(LowRankLDLT(L, T[1, 1e-4], check)) ≈ T(1e-4)
    end
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
