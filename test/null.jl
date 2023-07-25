module TestNull

using Test
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MultivariateMoments as MM

b(x) = MB.MonomialBasis(x)

include("utils.jl")

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
    @testset "$D" for D in [MM.AnyDependence, MM.StaircaseDependence]
        @testset "$name" for (solver, name) in [
            (MM.ShiftNullspace{D}(), "shift"),
            #(MM.Echelon{D}(fallback = false), "echelon no fallback"),
            #(MM.Echelon{D}(fallback = true), "echelon fallback"),
        ]
            @test isnothing(MM.solve(null, solver))
        end
    end
end

function _test_univariate_infinity(x, y)
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
        @testset "$D" for D in [MM.AnyDependence, MM.StaircaseDependence]
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
