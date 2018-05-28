@testset "Measure" begin
    Mod.@polyvar x y
    @test_throws ArgumentError Measure([1, 2], [x, x*y, y])
    @test_throws ArgumentError Measure([1, 2, 3, 4], [x, x*y, y])
    m = Measure([1, 0, 2, 3], [x^2*y^2, y*x^2, x*y*x^2, x*y^2])
    @test m.a == [2, 1, 0, 3]
    @test (-m).a == [-2, -1, 0, -3]
    @test (2 * m).a == [4, 2, 0, 6]
    @test (m * 3).a == [6, 3, 0, 9]
    #@test_throws ArgumentError Measure([1], [x]) + Measure([1], [y])
end
