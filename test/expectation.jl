@testset "Expectation" begin
    Mod.@polyvar x[1:3]
    p = x[3] - 2x[1]*x[2]^2 + 3x[3]*x[1] - 5x[1]^3
    v = (1,2,3)
    m = dirac(monomials(p), x=>v)
    @test MultivariateMoments.expectation(m, p) == p(x=>v) == MultivariateMoments.expectation(p, m)
    @test_throws ErrorException dot(x[1] * x[2] * x[3], m)
    @test dot(0.5 * x[1] * x[2]^2, m) == 2.0
    @test dot(m, x[1] * x[3]) == 3

    @testset "Expectation - MB" begin
        Mod.@polyvar x[1:2]
        p = x[2] - 2x[1]*x[2]^2 + 3x[2]*x[1] - 5x[1]^3
        m = measure([1, 0, 2, 3, 0, 4, 0, 2, 5, 0], maxdegree_basis(ChebyshevBasis, x, 3))
        @test MultivariateMoments.expectation(m, p) == MultivariateMoments.expectation(p, m) == 4.25
        @test_throws ErrorException dot((x[1] * x[2])^2, m)
        @test dot(0.5 * x[1] * x[2]^2, m) == 1.0
        @test dot(m, x[1] * x[2]) == 4 
    end

end
