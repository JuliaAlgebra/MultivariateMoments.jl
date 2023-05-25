@testset "Atomic measures" begin
    Mod.@polyvar x y
    δ1 = WeightedDiracMeasure([1, 0], 2.0)
    δ2 = WeightedDiracMeasure([1 / 2, 1 / 2], 3.0)
    η = AtomicMeasure([x, y], [δ1, δ2])
    @test η ≈ AtomicMeasure([x, y], [δ2, δ1])
    @test !(
        η ≈ AtomicMeasure(
            [x, y],
            [WeightedDiracMeasure([1, 0], 3.0), WeightedDiracMeasure([1 / 2, 1 / 2], 2.0)],
        )
    )
    @test string(η) ==
          "Atomic measure on the variables x, y with 2 atoms:\n at [1.0, 0.0] with weight 2.0\n at [0.5, 0.5] with weight 3.0"
    p = x^2 + x + y + x * y^2 + 1
    @test dot(η, p) == 2 * p(x => 1, y => 0) + 3 * p(x => 1 / 2, y => 1 / 2)
    @test dot(η, p) == dot(p, η)
    μ = measure(η, monomials(p))
    @test dot(μ, p) == 2 * p(x => 1, y => 0) + 3 * p(x => 1 / 2, y => 1 / 2)
end
