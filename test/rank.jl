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

struct HardcodedRanks <: RankCheck
    r::Vector{Int}
end

MultivariateMoments.rank_from_singular_values(σ, r::HardcodedRanks) = r.r[length(σ)]

@testset "Decreasing rank" begin
    r = RankDependence(zeros(4, 4), HardcodedRanks([1, 0, 1, 2]))
    @test !is_dependent!(r, 1)
    @test is_dependent!(r, 2)
    @test is_dependent!(r, 3)
    @test is_dependent!(r, 4)
end
