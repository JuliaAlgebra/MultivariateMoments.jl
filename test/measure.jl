@testset "Measure" begin
    Mod.@polyvar x y
    @test_throws ArgumentError measure([1, 2], [x, x*y, y])
    @test_throws ArgumentError measure([1, 2, 3, 4], [x, x*y, y])
    μ = measure([1, 0, 2, 3], [x^2*y^2, y*x^2, x*y*x^2, x*y^2])
    @test monomials(μ) == [x^3*y, x^2*y^2, x^2*y, x*y^2]
    @test monomial.(moments(μ)) == [x^3*y, x^2*y^2, x^2*y, x*y^2]
    @test MultivariateMoments.moment_value.(moments(μ)) == [2, 1, 0, 3]
    @test all(nvariables.(moments(μ)) .== 2)
    @test degree.(moments(μ)) == [4, 4, 3, 3]
    @test μ.a == [2, 1, 0, 3]
    @test (-μ).a == [-2, -1, 0, -3]
    @test (2 * μ).a == [4, 2, 0, 6]
    @test (μ * 3).a == [6, 3, 0, 9]
    #@test_throws ArgumentError measure([1], [x]) + measure([1], [y])
end
