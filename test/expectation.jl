@testset "Expectation" begin
    Mod.@polyvar x[1:3]
    p = x[3] - 2x[1] * x[2]^2 + 3x[3] * x[1] - 5x[1]^3
    v = (1, 2, 3)
    m = dirac(monomials(p), x => v)
    @test (@inferred MultivariateMoments.expectation(m, p)) ==
          p(x => v) ==
          (@inferred MultivariateMoments.expectation(p, m))
    @test_throws ErrorException dot(x[1] * x[2] * x[3], m)
    @test (@inferred dot(0.5 * x[1] * x[2]^2, m)) == 2.0
    @test (@inferred dot(m, x[1] * x[3])) == 3
end
