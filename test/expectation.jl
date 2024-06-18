@testset "Expectation" begin
    Mod.@polyvar x[1:3]
    p = x[3] - 2x[1] * x[2]^2 + 3x[3] * x[1] - 5x[1]^3
    v = (1, 2, 3)
    m = dirac(monomials(p), x => v)
    a = MB.algebra_element(p)
    expected = p(x => v)
    @test (@inferred MultivariateMoments.expectation(m, p)) == expected
    @test (@inferred MultivariateMoments.expectation(m, a)) == expected
    @test (@inferred MultivariateMoments.expectation(p, m)) == expected
    @test (@inferred MultivariateMoments.expectation(a, m)) == expected
    @test (@inferred dot(m, p)) == expected
    @test (@inferred dot(m, a)) == expected
    @test (@inferred dot(p, m)) == expected
    @test (@inferred dot(a, m)) == expected
    @test_throws ErrorException dot(x[1] * x[2] * x[3], m)
    @test (@inferred dot(0.5 * x[1] * x[2]^2, m)) == 2.0
    a = MB.algebra_element([0.5], MB.SubBasis{MB.Monomial}([x[1] * x[2]^2]))
    @test (@inferred dot(a, m)) == 2.0
    a = MB.algebra_element([0.5], MB.SubBasis{MB.ScaledMonomial}([x[1] * x[2]^2]))
    @test (@inferred dot(a, m)) == 2√3
    @test (@inferred dot(m, x[1] * x[3])) == 3
end
