using Random, Test
using MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
using SemialgebraicSets
using MultivariateMoments

struct DummySolver <: SemialgebraicSets.AbstractAlgebraicSolver end
function SemialgebraicSets.solve(
    ::SemialgebraicSets.AlgebraicSet,
    ::DummySolver,
)
    error("Dummy solver")
end

# [HL05] Henrion, D. & Lasserre, J-B.
# Detecting Global Optimality and Extracting Solutions of GloptiPoly
# 2005
#
# [LPJ20] Legat, Benoît, Pablo Parrilo, and Raphaël Jungers.
# "Certifying unstability of switched systems using sum of squares programming."
# SIAM Journal on Control and Optimization 58.4 (2020): 2616-2638.

Mod.@polyvar x

function testelements(X, Y, atol)
    @test length(X) == length(Y)
    for y in Y
        @test any(x -> isapprox(x, y, atol=atol), X)
    end
end

function _atoms(atoms, rank_check, solver)
    Mod.@polyvar x[1:length(first(atoms))]
    η = AtomicMeasure(x, WeightedDiracMeasure.(atoms, ones(length(atoms))))
    monos = monomials(x, 0:(length(atoms) + 2))
    μ = measure(η, monos)
    ν = moment_matrix(μ, monomials(x, 0:(div(length(atoms), 2) + 1)))
    atoms = atomic_measure(ν, rank_check, solver)
    @test atoms !== nothing
    @test atoms ≈ η
end

function atoms_1(rank_check, solver)
    atoms = [
        [1, 2],
    ]
    _atoms(atoms, rank_check, solver)
    return
end

function atoms_2(rank_check, solver)
    atoms = [
        [1, 1],
        [1, -1],
    ]
    _atoms(atoms, rank_check, solver)
    return
end
"""
    hl05_2_3()

[HL05, Section 2.3]
"""
function hl05_2_3(rank_check, lrc, solver, perturb::Bool=true)
    Mod.@polyvar x y
    η = AtomicMeasure([x, y], WeightedDiracMeasure.([[1, 2], [2, 2], [2, 3]], [0.4132, 0.3391, 0.2477]))
    monos = [x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1]
    μ = measure(η, monos)
    ν = moment_matrix(μ, [1, x, y, x^2, x*y, y^2])
    atoms = atomic_measure(ν, rank_check, lrc, Echelon(), solver)
    @test atoms !== nothing
    @test atoms ≈ η
    if perturb # the shift `1e-14` is too small compared to the noise of `1e-6`. We want high noise so that the default rtol of `Base.rtoldefault` does not work so that it tests that `rtol` is passed around.
        Random.seed!(0)
        ν2 = MomentMatrix(SymMatrix(ν.Q.Q + rand(length(ν.Q.Q)) * 1e-6, ν.Q.n), ν.basis)
        @test_throws ErrorException atomic_measure(ν2, rank_check, lrc, Echelon(), solver, weight_solver = MomentVectorWeightSolver())
        for weight_solver in [MomentMatrixWeightSolver(), MomentVectorWeightSolver(rtol=1e-5), MomentVectorWeightSolver(atol=1e-5)]
            atoms = atomic_measure(ν2, rank_check, lrc, Echelon(), solver; weight_solver)
            @test atoms !== nothing
            @test atoms ≈ η rtol=1e-4
        end
    end
end

"""
    hl05_3_3_1()

[HL05, Section 3.3.1]
"""
function hl05_3_3_1()
    Mod.@polyvar x y z
    U = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         2 0 0 0 0 0; # z^2 = 2z
         0 0 0 1 0 0;
         0 2 0 0 0 0; # y^2 = 2y
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         0 0 2 0 0 0] # x^2 = 2x
    # β will be [z, y, x, y*z, x*z, x*y]
    x = monomial_vector([z, y, x, z^2, y*z, y^2, x*z, x*y, x^2])
    # x*β contains x*y*z, x^2*z, x^2*y which are not present so it show fail
    V = MultivariateMoments.build_system(U', MB.MonomialBasis(x), sqrt(eps(Float64)))
    @test is_zero_dimensional(V)
    testelements(V, [[2.0, 2.0, 2.0], [2.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 0.0]], 1e-11)
end

"""
    hl05_4()

[HL05, Section 4]
"""
function hl05_4(rank_check, lrc)
    Mod.@polyvar x y
    s3 = sqrt(1/3)
    η = AtomicMeasure([x, y], WeightedDiracMeasure.([[-s3, -s3], [-s3, s3], [s3, -s3], [s3, s3]], [0.25, 0.25, 0.25, 0.25]))
    μ = measure([1/9,     0,     1/9,     0, 1/9,   0,     0,     0,   0, 1/3,   0, 1/3, 0, 0, 1],
                [x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1])
    ν = moment_matrix(μ, [1, x, y, x^2, x*y, y^2])
    atoms = atomic_measure(ν, rank_check, lrc)
    @test atoms !== nothing
    if lrc isa LowRankChol
        @test atoms ≈ η
    end
    @test measure(atoms, μ.x).a ≈ μ.a rtol=1e-3
end

"""
    lpj20_3_8_0()

[LPJ20, Example 3.8]
"""
function lpj20_3_8_0(rank_check, lrc, ok::Bool=true)
    Mod.@polyvar x[1:2]
    η = AtomicMeasure(x, [WeightedDiracMeasure([1., 0.], 2.)])
    μ = measure([1e-6, 0.0, 2],
                monomials(x, 2))
    ν = moment_matrix(μ, monomials(x, 1))
    atoms = atomic_measure(ν, rank_check, lrc)
    if ok
        @test atoms !== nothing
        @test atoms ≈ η
    else
        @test atoms === nothing
    end
end

"""
    lpj20_3_8()

[LPJ20, Example 3.8]
"""
function lpj20_3_8(rank_check, solver)
    Mod.@polyvar x[1:2]
    η = AtomicMeasure(x, [WeightedDiracMeasure([1., 0.], 2.53267)])
    ν = moment_matrix(
        [ 2.53267 -0.0 -5.36283e-19
         -0.0 -5.36283e-19 -0.0
         -5.36283e-19 -0.0 7.44133e-6][3:-1:1, 3:-1:1],
        monomials(x, 2),
    )
    atoms = atomic_measure(ν, rank_check, Echelon(), solver)
    @test atoms !== nothing
    @test atoms ≈ η
end

"""
    lpj20_3_9()

[LPJ20, Example 3.9]
"""
function lpj20_3_9(rank_check, sol=1)
    Mod.@polyvar x[1:2]
    # η1 is slightly closer to ν than η2
    ν = moment_matrix([
        14.8107 3.37426 0.769196 0.175447 0.0400656 0.00917416 0.00211654
        3.37426 0.769196 0.175447 0.0400656 0.00917416 0.00211654 0.000498924
        0.769196 0.175447 0.0400656 0.00917416 0.00211654 0.000498924 0.000125543
        0.175447 0.0400656 0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5
        0.0400656 0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5
        0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5 1.07817e-5
        0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5 1.07817e-5 1.10658e-5
    ][7:-1:1, 7:-1:1], monomials(x, 6))
    atoms = atomic_measure(ν, rank_check, weight_solver=MomentVectorWeightSolver())
    if sol == 0
        @test atoms === nothing
    else
        @test atoms !== nothing
        if sol == 1
            η = AtomicMeasure(x, [WeightedDiracMeasure([1., 0.2278], 14.8106)])
        else
            η = AtomicMeasure(x, [WeightedDiracMeasure([1., 11.15240711893039],   6.038060135782876e-19),
                                  WeightedDiracMeasure([1., 0.991590438508714],   8.48980575400066e-6),
                                  WeightedDiracMeasure([1., 0.5482585022039119],  0.0011737977570447154),
                                  WeightedDiracMeasure([1., -0.2808541188399676], 9.72284626287903e-5),
                                  WeightedDiracMeasure([1., 0.22799858558381308], 14.781765717198045),
                                  WeightedDiracMeasure([1., 0.12343866254667915], 0.02765476684449834)])
        end
        @test atoms ≈ η rtol=1e-4
    end
end

"""
    jcg14_6_1()

[LPJ20] applied to [JCG14, Example 6.1].
"""
function jcg14_6_1(rank_check, ok::Bool=true)
    Mod.@polyvar x[1:4]
    η = AtomicMeasure(x, [WeightedDiracMeasure([1.0, 4.78736282579504, 1.24375738760842, -1.4231836829978], 0.039112791390926646)])
    ν = moment_matrix(
        [ 0.0397951 0.187094 0.0489553 -0.0551816
          0.187094 0.896353 0.232962 -0.265564
          0.0489553 0.232962 0.0614682 -0.0676226
         -0.0551816 -0.265564 -0.0676226 0.0837186][4:-1:1, 4:-1:1],
        monomials(x, 1),
    )
    atoms = atomic_measure(ν, rank_check)
    if ok
        @test atoms !== nothing
        @test atoms ≈ η
    else
        @test atoms === nothing
    end
end

function large_norm(rank_check)
    # If the norm of `M` is given to `rref!` instead of `√||M||`, `atomic_measure` will error.
    Mod.@polyvar x[1:2]
    Q = [
        586.8034549325414 -800.152792847183
       -800.152792847183  2749.376669556701
    ]
    ν = moment_matrix(Q, monomials(x, 1))
    @test nothing === atomic_measure(ν, rank_check)
end

function test_extract()
    default_solver = SemialgebraicSets.default_algebraic_solver([1.0x - 1.0x])
    for solver in [SVDCholesky(), ShiftCholesky(1e-15), FlatExtension(), FlatExtension(IterativeDiagonalization()), ShiftNullspace()]
        atoms_1(1e-10, solver)
        atoms_2(1e-10, solver)
    end
    for lrc in (SVDCholesky(), ShiftCholesky(1e-14))
        perturb = !(lrc isa ShiftCholesky) # the shift `1e-14` is too small compared to the noise of `1e-6`. We want high noise so that the default rtol of `Base.rtoldefault` does not work so that it tests that `rtol` is passed around.
        hl05_2_3(1e-4, lrc, default_solver, perturb)
        @test_throws ErrorException("Dummy solver") hl05_2_3(1e-4, lrc, DummySolver())
    end
    hl05_3_3_1()
    # Fails on 32-bits in CI
    if Sys.WORD_SIZE != 32
        for lrc in (SVDCholesky(), ShiftCholesky(1e-16))
            hl05_4(1e-16, lrc)
        end
    end
    # All singular values will be at least 1e-6 > 1e-12 it won't eliminate any row
    lpj20_3_8_0(1e-12, ShiftCholesky(1e-6), false)
    # The following tests that the method does not error if ranktol eliminates everything
    # In particular, this tests that the function equation(i) do not call sum when r equal to 0
    # this that throws an ArgumentError as details in src/extract.jl
    lpj20_3_8_0(1.0, SVDCholesky(), false)
    lpj20_3_8_0(LeadingRelativeRankTol(1e-5), SVDCholesky())
    lpj20_3_8_0(AbsoluteRankTol(1e-5), SVDCholesky())
    lpj20_3_8_0(DifferentialRankTol(1e-2), SVDCholesky())
    lpj20_3_8_0(LargestDifferentialRank(), SVDCholesky())
    @test_throws ErrorException("Dummy solver") lpj20_3_8(1e-5, DummySolver())
    lpj20_3_8(1e-5, default_solver)
    for ranktol in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
        # With 1e-3 and 1e-4, the rank is detected to be 1
        # With 1e-5, the rank is detected to be 2
        # With 1e-6, the rank is detected to be 3
        # With 1e-7, the rank is detected to be 5
        lpj20_3_9(ranktol)
    end
    for fixed_rank in 1:5
        lpj20_3_9(FixedRank(fixed_rank))
    end
    # With 1e-8, the rank is detected to be 6
    lpj20_3_9(1e-8, 2)
    lpj20_3_9(FixedRank(6), 2)
    # With 1e-9, the rank is detected to be 7
    lpj20_3_9(1e-9, 0)
    lpj20_3_9(FixedRank(7), 0)
    jcg14_6_1(6e-3)
    jcg14_6_1(8e-4, false)
    large_norm(1e-2)
end

@testset "Atom extraction" begin
    test_extract()
end
