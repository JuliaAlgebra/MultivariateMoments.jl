using Test, LinearAlgebra, MultivariateMoments
import PositiveFactorizations
import DynamicPolynomials as DP
import MultivariatePolynomials as MP

@testset "PositiveFactorizations extension" begin
    @test Base.get_extension(
        MultivariateMoments,
        :MultivariateMomentsPositiveFactorizationsExt,
    ) !== nothing
    algo = AbsCholeskyLDLT()
    for T in (Float32, Float64, BigFloat)
        tol = sqrt(eps(T))
        for (A, expected_rank) in (
            (T[1 -1; -1 1], 1),
            (T[0 0 0; 0 4 2; 0 2 1], 1),
            (T[1 0 0; 0 9 0; 0 0 4], 3),
            (zeros(T, 2, 2), 0),
            (zeros(T, 0, 0), 0),
        )
            original = copy(A)
            check = LeadingRelativeRankTol(tol)
            F = low_rank_ldlt(Symmetric(A), algo, check)
            r = size(F.L, 2)
            @test F.L * Diagonal(F.singular_values[1:r]) * F.L' ≈ A
            @test r == expected_rank
            @test issorted(F.singular_values; rev = true)
            @test F.rank_check === check
            @test all(isfinite, F.L)
            @test eltype(F.L) == T
            @test A == original
        end
    end
    F = low_rank_ldlt([0 0; 0 4], algo, 1e-10)
    @test F.singular_values == [4.0, 0.0]
    @test F.L * Diagonal(F.singular_values[1:1]) * F.L' ≈ [0 0; 0 4]
    # The negative pivot is made positive instead of requiring a shift.
    F = low_rank_ldlt([-4.0 2.0; 2.0 -1.0], algo, 1e-10)
    @test F.L * Diagonal(F.singular_values[1:1]) * F.L' ≈ [4 -2; -2 1]
    @test size(F.L) == (2, 1)
    F = low_rank_ldlt([1.0 0.0; 0.0 1e-8], algo, 1e-6)
    @test size(F.L) == (2, 1)
    @test accuracy(F) ≈ 1e-8
    @test size(low_rank_ldlt(Matrix{Float64}(I, 2, 2), algo, FixedRank(1)).L) ==
          (2, 1)
    DP.@polyvar x
    η = AtomicMeasure([x], [WeightedDiracMeasure([2.0], 3.0)])
    ν = moment_matrix(
        moment_vector(η, MP.monomials(x, 0:4)),
        MP.monomials(x, 0:2),
    )
    @test atomic_measure(ν, LeadingRelativeRankTol(1e-10), algo) ≈ η
    @test_throws DimensionMismatch low_rank_ldlt(zeros(2, 3), algo, 1e-10)
end
