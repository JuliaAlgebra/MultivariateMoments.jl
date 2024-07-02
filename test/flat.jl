using Test, MultivariateMoments

@testset "ZeroDimensionalVariety" begin
    V = ZeroDimensionalVariety([[1], [2]])
    expected = "ZeroDimensionalVariety with elements:\n[[1], [2]]"
    @test sprint(show, V) == expected
    @test sprint(show, MIME"text/plain"(), V) == expected
end
