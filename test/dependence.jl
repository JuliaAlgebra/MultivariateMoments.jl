module TestDependence

using Test
import RecipesBase as RB
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MultivariateMoments as MM

b(x) = MB.MonomialBasis(x)

struct FixedDependence
    dependent::Vector{Int}
end

MM.is_dependent!(d::FixedDependence, i) = i in d.dependent
function MM.column_compression!(::FixedDependence, _) end

function test_degree_error(x, y, z)
    M = typeof(x * y * z)
    m = b(M[])
    a = MM.AnyDependence(_ -> true, m)
    @test sprint(show, a) == "AnyDependence for an empty basis\n"
    @test isempty(a)
    s = MM.StaircaseDependence(_ -> true, m)
    @test sprint(show, s) == "StaircaseDependence for an empty basis\n"
    @test isempty(s)
#    for (f, name) in [(MP.mindegree, "min"), (MP.maxdegree, "max")]
#        err = ErrorException(
#            "Cannot compute `$(name)degree` as all bases are empty",
#        )
#        @test_throws err f(a)
#        @test_throws err f(a, x)
#        @test_throws err f(s)
#        @test_throws err f(s, x)
#    end
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
    basis = b([1, x^2])
    d = MM.StaircaseDependence(FixedDependence([2]), basis)
    @test d.basis.monomials == [x^0, x, x^2]
    @test d.dependence == [
        MM.LinearDependence(false, true, false),
        MM.LinearDependence(false, false, false),
        MM.LinearDependence(true, true, false),
    ]
    @test d.position == [MM.STANDARD, MM.STANDARD, MM.CORNER]
    d = MM.StaircaseDependence(FixedDependence(Int[]), b([1, x]))
    @test d.basis.monomials == [x^0, x, x^2]
    @test d.dependence == [
        MM.LinearDependence(false, true, false),
        MM.LinearDependence(false, true, false),
        MM.LinearDependence(false, false, false),
    ]
    @test d.position == [MM.STANDARD, MM.STANDARD, MM.STANDARD]
end

function test_recipe(x, y, z)
    a = [x^0 * y^0]
    A = ([0], [0])
    c = [x, y^2]
    C = ([1, 0], [0, 2])
    ac, Ia, Ic = MB.merge_bases(b(a), b(c))
    d = [x * z^1, x * y^2]
    D = ([1, 1], [0, 2], [1, 0])
    e = [x * y^2, x * y^3, x * z^4]
    E = ([1, 1, 1], [2, 3, 0], [0, 0, 4])
    ae, _, _ = MB.merge_bases(b(a), b(e))
    cd, _, _ = MB.merge_bases(b(c), b(d))
    aecd, _, Idep = MB.merge_bases(ae, cd)
    _test_recipe(
        MM.AnyDependence(i -> !iszero(Ic), ac),
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
    _test_recipe(
        MM.StaircaseDependence(
            FixedDependence(findall(!iszero, Idep)),
            aecd,
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
    return
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
