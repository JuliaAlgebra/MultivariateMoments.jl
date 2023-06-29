using Test

@testset "MomentMatrix" begin
    Mod.@polyvar x y
    μ = measure([1], [x])
    for mono in [x^0, x, y]
        err = ArgumentError("`$μ` does not have the moment `$(mono^2)`")
        @test_throws err moment_matrix(μ, [mono])
    end
    μ = measure(1:3, [x^2, x * y, y^2])
    X = [x, y]
    ν1 = moment_matrix(μ, X)
    @test sprint(show, ν1) == """
MomentMatrix with row/column basis:
 MonomialBasis([y, x])
And entries in a 2×2 SymMatrix{Int64}:
 3  2
 2  1"""
    ν2 = MomentMatrix{Int}((i, j) -> i + j - 1, X)
    for ν in (ν1, ν2)
        @test ν.Q[1:4] == [3, 2, 2, 1]
        @test ν.Q[1, 1] == 3
        @test ν.Q[1, 2] == 2
        @test ν.Q[2, 1] == 2
        @test ν.Q[2, 2] == 1
        @test variables(ν)[1] == x
        @test variables(ν)[2] == y
        @test nvariables(ν) == 2
    end
    @test_throws ArgumentError moment_matrix(measure([1], [x]), [y])
    block_ν = BlockDiagonalMomentMatrix([ν1, ν2])
    @test block_ν isa BlockDiagonalMomentMatrix{Int,typeof(ν1.basis)}
    @test sprint(show, block_ν) == """
BlockDiagonalMomentMatrix with 2 blocks:
[1] Block with row/column basis:
     MonomialBasis([y, x])
    And entries in a 2×2 SymMatrix{Int64}:
     3  2
     2  1
[2] Block with row/column basis:
     MonomialBasis([y, x])
    And entries in a 2×2 SymMatrix{Int64}:
     3  2
     2  1"""
end
