@testset "low_rank_ldlt" begin
    v = [1, -1]
    M = v * v'
    ldlt = low_rank_ldlt(M, SVDLDLT(), 1e-10)
    @test ldlt.L ≈ -normalize(v)
    ldlt = low_rank_ldlt(M, ShiftCholeskyLDLT(1e-10), 1e-9)
    @test ldlt.L ≈ v
end
