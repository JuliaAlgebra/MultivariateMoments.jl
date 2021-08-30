@testset "Measure" begin
    Mod.@polyvar x y
    @test_throws ArgumentError measure([1, 2], [x, x*y, y])
    @test_throws ArgumentError measure([1, 2, 3, 4], [x, x*y, y])
    μ = measure([1, 0, 2, 3], [x^2*y^2, y*x^2, x*y*x^2, x*y^2])
    @test base_functions(μ) == MonomialBasis(monovec([x^3*y, x^2*y^2, x^2*y, x*y^2]))
    @test base_function.(moments(μ)) == [x^3*y, x^2*y^2, x^2*y, x*y^2]
    @test MultivariateMoments.moment_value.(moments(μ)) == [2, 1, 0, 3]
    @test all(nvariables.(moments(μ)) .== 2)
    @test degree.(moments(μ)) == [4, 4, 3, 3]
    @test μ.values == [2, 1, 0, 3]
    @test (-μ).values == [-2, -1, 0, -3]
    @test (2 * μ).values == [4, 2, 0, 6]
    @test (μ * 3).values == [6, 3, 0, 9]
    #@test_throws ArgumentError measure([1], [x]) + measure([1], [y])

    @testset "Measure MB" begin
        Mod.@polyvar x y
        @test_throws AssertionError measure([1, 0, 2, 3, 0], maxdegree_basis(ChebyshevBasis, [x, y], 2))
        @test_throws AssertionError measure([1, 0, 2, 3, 0, 4], maxdegree_basis(ChebyshevBasis, [x, y], 3))
        μ = measure([1, 0, 2, 3, 0, 4], maxdegree_basis(ChebyshevBasis, [x, y], 2))
        @test base_functions(μ) ==  maxdegree_basis(ChebyshevBasis, [x, y], 2)
        @test base_function.(moments(μ)) == [2.0x^2 - 1.0, x*y, 2.0y^2 - 1.0, x, y, 1.0]
        @test MultivariateMoments.moment_value.(moments(μ)) == [1, 0, 2, 3, 0, 4]
        @test all(nvariables.(moments(μ)) .== 2)
        @test degree.(moments(μ)) == [2, 2, 2, 1, 1, 0]
        @test μ.values == [1, 0, 2, 3, 0, 4]
        @test (-μ).values == -[1, 0, 2, 3, 0, 4]
        @test (2 * μ).values == [1, 0, 2, 3, 0, 4].*2
        @test (μ * 3).values == [1, 0, 2, 3, 0, 4].*3

    end
end
