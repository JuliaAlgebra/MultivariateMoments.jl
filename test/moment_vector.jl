import StarAlgebras as SA

@testset "MomentVector" begin
    Mod.@polyvar x y
    @test_throws DimensionMismatch moment_vector([1, 2], [x, x * y, y])
    @test_throws DimensionMismatch moment_vector([1, 2, 3, 4], [x, x * y, y])
    μ = moment_vector([1, 0, 2, 3], [x^2 * y^2, y * x^2, x * y * x^2, x * y^2])
    @test MP.mindegree(μ) == 3
    @test MP.maxdegree(μ) == 4
    @test MP.extdegree(μ) == (3, 4)
    @test SA.basis(μ).monomials ==
          monomial_vector([x^3 * y, x^2 * y^2, x^2 * y, x * y^2])
    @test map(m -> m.polynomial.monomial, moments(μ)) ==
          monomial_vector([x^3 * y, x^2 * y^2, x^2 * y, x * y^2])
    @test MultivariateMoments.moment_value.(moments(μ)) == reverse([2, 1, 0, 3])
    @test all(m -> variables(m) == variables(x * y), moments(μ))
    @test all(nvariables.(moments(μ)) .== 2)
    @test μ.values == reverse([2, 1, 0, 3])
    @test (-μ).values == reverse([-2, -1, 0, -3])
    @test (2 * μ).values == reverse([4, 2, 0, 6])
    @test (μ * 3).values == reverse([6, 3, 0, 9])
    μ = moment_vector([1, 1], [x, x])
    @test MultivariateMoments.moment_value.(moments(μ)) == [1]
    err = ErrorException(
        "The monomial `x` occurs twice with different values: `1` and `2`",
    )
    @test_throws err moment_vector([1, 2], [x, x])
    #@test_throws ArgumentError moment_vector([1], [x]) + moment_vector([1], [y])
end
