module TestNull

using Test
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MultivariateMoments as MM

b(x) = MB.MonomialBasis(x)

include("utils.jl")

function test_partial_commutative_fix(x, y)
    matrix = Float64[
        1 0 0 0 # 1
        0 1 0 0 # y
        0 0 1 0 # x
        0 0 0 1 # y^2
        0 0 1 0 # x * y
        1 0 0 0 # x^2
    ]
    monos = MP.monomials((x, y), 0:2)
    for basis in [MB.MonomialBasis(monos), MB.OrthonormalCoefficientsBasis(monos)]
        null = MM.MacaulayNullspace(matrix, basis, 1e-8)
        D = MM.StaircaseDependence
        solver = MM.StaircaseSolver{Float64}(max_partial_iterations = 1)
        shift_solver = MM.ShiftNullspace{D}(solver)
        sol = MM.solve(null, shift_solver)
        testelements(sol, [[1, 1], [-1, 1]], 1e-8)
    end
    return
end

function test_dependent_border(x, y)
    matrix = Float64[
        1 0 0 0 0 # 1
        0 1 0 0 0 # y
        0 0 1 0 0 # x
        0 0 0 1 0 # y^2
        -1 0 0 0 0 # x * y = -1
        0 0 0 0 1 # x^2
        0 1 0 0 0 # y^3 = y
        0 0 1 0 0 # x * y^2 = x
        0 1 0 0 0 # x^2 * y = y
        0 0 1 0 0 # x^3 = x
    ]
    basis = MB.MonomialBasis(MP.monomials((x, y), 0:3))
    null = MM.MacaulayNullspace(matrix, basis, 1e-8)
    @testset "$D" for D in [MM.LinearDependence, MM.StaircaseDependence]
        @testset "$name" for (solver, name) in [
            (MM.ShiftNullspace{D}(), "shift"),
            (MM.Echelon{D}(fallback = false), "echelon no fallback"),
            (MM.Echelon{D}(fallback = true), "echelon fallback"),
        ]
            sol = MM.solve(null, solver)
            testelements(sol, [[1, -1], [-1, 1]], 1e-8)
        end
    end
end

function test_independent_border(x, y)
    matrix = Float64[
        1 0 0 0 0 0 # 1
        0 1 0 0 0 0 # y
        0 0 1 0 0 0 # x
        0 0 0 1 0 0 # y^2
        -1 0 0 0 0 0 # x * y = -1
        0 0 0 0 1 0 # x^2
        0 1 0 0 0 0 # y^3 = y
        0 0 1 0 0 0 # x * y^2 = x
        0 0 0 0 0 1 # x^2 * y: indepedent border
        0 0 1 0 0 0 # x^3 = x
    ]
    basis = MB.MonomialBasis(MP.monomials((x, y), 0:3))
    null = MM.MacaulayNullspace(matrix, basis, 1e-8)
    @testset "$D" for D in [MM.LinearDependence, MM.StaircaseDependence]
        @testset "$name" for (solver, name) in [
            (MM.ShiftNullspace{D}(), "shift"),
            (MM.Echelon{D}(fallback = false), "echelon no fallback"),
            (MM.Echelon{D}(fallback = true), "echelon fallback"),
        ]
            sol = MM.solve(null, solver)
            testelements(sol, [[1, -1], [-1, 1]], 1e-8)
        end
    end
end

function test_dreesen1(x, y)
    matrix = Float64[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1
        0 0 1 0
        1 0 0 0
    ]
    basis = MB.MonomialBasis([1, y, x, y^2, x * y, x^2])
    null = MM.MacaulayNullspace(matrix, basis, 1e-8)
    @testset "$D" for D in [MM.LinearDependence, MM.StaircaseDependence]
        @testset "$name" for (solver, name) in [
            (MM.ShiftNullspace{D}(), "shift"),
            #(MM.Echelon{D}(fallback = false), "echelon no fallback"),
            #(MM.Echelon{D}(fallback = true), "echelon fallback"),
        ]
            @test isnothing(MM.solve(null, solver))
        end
    end
end

function test_univariate_infinity(x, y)
    Z = Float64[
        1 0
        2 0
        0 1
    ]
    for (null, sol) in [
        (MM.MacaulayNullspace(Z, b([1, x, x^2]), 1e-8), 2),
        (MM.MacaulayNullspace(Z[1:2, 1:1], b([1, x]), 1e-8), 2),
        (MM.MacaulayNullspace(zeros(1, 0), b([x]), 1e-8), 0),
    ]
        @testset "$D" for D in [MM.LinearDependence, MM.StaircaseDependence]
            @testset "$name" for (solver, name) in [
                (MM.ShiftNullspace{D}(), "shift"),
                (MM.Echelon{D}(fallback = false), "echelon no fallback"),
                (MM.Echelon{D}(fallback = true), "echelon fallback"),
            ]
                if name == "echelon no fallback" && size(null.matrix) == (1, 0)
                    # FIXME The system is `x = 0` hence it thing that it is homogeneous and it finds no sol
                    #       The homogeneous check should be better... It should probably be explicitly stated upfront
                    continue
                end
                @test collect(MM.solve(null, solver)) == [[sol]]
            end
        end
    end
end

function runtests(args...)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)(args...)
            end
        end
    end
end

end

using Test

import DynamicPolynomials
@testset "DynamicPolynomials" begin
    DynamicPolynomials.@polyvar x y
    TestNull.runtests(x, y)
end

import TypedPolynomials
@testset "TypedPolynomials" begin
    TypedPolynomials.@polyvar x y
    TestNull.runtests(x, y)
end
