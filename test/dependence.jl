module TestDependence

using Test
import RecipesBase as RB
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MultivariateMoments as MM

b(x) = MB.MonomialBasis(x)

function test_degree_error(x, y, z)
    M = typeof(x * y * z)
    m = b(M[])
    a = MM.AnyDependence(m, m)
    @test sprint(show, a) == "AnyDependence\n"
    @test isempty(a)
    s = MM.StaircaseDependence(m, m, m, m, m)
    @test sprint(show, s) == "StaircaseDependence\n"
    @test isempty(s)
    for (f, name) in [(MP.mindegree, "min"), (MP.maxdegree, "max")]
        err = ErrorException(
            "Cannot compute `$(name)degree` as all bases are empty",
        )
        @test_throws err f(a)
        @test_throws err f(a, x)
        @test_throws err f(s)
        @test_throws err f(s, x)
    end
end

function _test_recipe(dep, ticks, args, names, indep)
    @test sprint(show, dep) isa String
    d = Dict{Symbol,Any}()
    r = RB.apply_recipe(d, dep)
    @test length(d) == 1 + length(ticks)
    @test d[:aspect_ratio] == :equal
    @test d[:xticks] == ticks[1]
    @test d[:yticks] == ticks[2]
    if length(ticks) >= 3
        @test d[:zticks] == ticks[3]
    end
    @test length(r) == length(names)
    for i in eachindex(names)
        @test r[i].args == args[i]
        @test r[i].plotattributes[:seriestype] == :scatter
        shape = indep[i] ? :circle : :rect
        @test r[i].plotattributes[:markershape] == shape
        @test r[i].plotattributes[:label] == names[i]
    end
end

function test_staircase(x, y, z)
    d = MM.StaircaseDependence(b([1, x^2])) do i
        return i > 1
    end
    @test d.trivial_standard.monomials == [x]
    @test d.standard.monomials == [x^0]
    @test d.corners.monomials == [x^2]
    @test isempty(d.dependent_border.monomials)
    @test isempty(d.independent_border.monomials)
end

function test_recipe(x, y, z)
    a = [x^0 * y^0]
    A = ([0], [0])
    c = [x, y^2]
    C = ([1, 0], [0, 2])
    d = [x * z^1, x * y^2]
    D = ([1, 1], [0, 2], [1, 0])
    e = [x * y^2, x * y^3, x * z^4]
    E = ([1, 1, 1], [2, 3, 0], [0, 0, 4])
    _test_recipe(
        MM.AnyDependence(b(a), b(c)),
        [0:1, 0:2],
        [A, C],
        ["Independent", "Dependent"],
        [true, false],
    )
    _test_recipe(
        MM.AnyDependence(b(d), b(e)),
        [1:1, 0:3, 0:4],
        [D, E],
        ["Independent", "Dependent"],
        [true, false],
    )
    return _test_recipe(
        MM.StaircaseDependence(
            b(a .* z^0),
            b([x^0 * y^0 * z]),
            b(c .* z^0),
            b(d),
            b(e),
        ),
        [0:1, 0:3, 0:4],
        [(A..., [0]), ([0], [0], [1]), (C..., [0, 0]), D, E],
        [
            "Trivial standard",
            "Standard",
            "Corners",
            "Dependent border",
            "Independent border",
        ],
        [true, true, false, false, true],
    )
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
    DynamicPolynomials.@polyvar x y z
    TestDependence.runtests(x, y, z)
end

import TypedPolynomials
@testset "TypedPolynomials" begin
    TypedPolynomials.@polyvar x y z
    TestDependence.runtests(x, y, z)
end
