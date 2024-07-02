using Test, MultivariateMoments

@testset "ZeroDimensionalVariety" begin
    V = ZeroDimensionalVariety([[1], [2]])
    @show sprint(show, V)
    @show sprint(show, MIME"text/plain"(), V)
end
