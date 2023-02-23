using Test

@testset "MomentMatrix" begin
    Mod.@polyvar x y
    μ = measure([1], [x])
    for mono in [x^0, x, y]
        err = ArgumentError("`$μ` does not have the moment `$(mono^2)`")
        @test_throws err moment_matrix(μ, [mono])
    end
    μ = measure(1:3, [x^2, x*y, y^2])
    X = [x, y]
    ν1 = moment_matrix(μ, X)
    ν2 = MomentMatrix{Int}((i, j) -> i + j - 1, X)
    for ν in (ν1, ν2)
        @test ν.Q[1:4] == [1, 2, 2, 3]
        @test ν.Q[1, 1] == 1
        @test ν.Q[1, 2] == 2
        @test ν.Q[2, 1] == 2
        @test ν.Q[2, 2] == 3
        @test variables(ν)[1] == x
        @test variables(ν)[2] == y
        @test nvariables(ν) == 2
    end
    @test_throws ArgumentError moment_matrix(measure([1], [x]), [y])
    sparse_ν = SparseMomentMatrix([ν1, ν2])
    @test sparse_ν isa SparseMomentMatrix{Int,typeof(ν1.basis)}
end
