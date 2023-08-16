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
    @test sprint(show, a) == "AnyDependence for an empty basis
"
    @test isempty(a)
    s = MM.StaircaseDependence(_ -> true, m)
    @test sprint(show, s) == "StaircaseDependence for an empty basis
"
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

function _test_recipe(dep, ticks, args, names, shapes)
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
    @testset "$i" for i in eachindex(names)
        @test r[i].args == args[i]
        @test r[i].plotattributes[:seriestype] == :scatter
        @test r[i].plotattributes[:markershape] == shapes[i]
        @test r[i].plotattributes[:label] == names[i]
    end
end

function test_staircase(x, y, z)
    basis = b([1, x^2])
    d = MM.StaircaseDependence(FixedDependence([2]), basis)
    @test d.basis.monomials == [x^0, x, x^2]
    @test d.dependence == [
        MM.INDEPENDENT,
        MM.TRIVIAL,
        MM.DEPENDENT,
    ]
    @test d.position == [MM.STANDARD, MM.STANDARD, MM.CORNER]
    d = MM.StaircaseDependence(FixedDependence(Int[]), b([1, x]))
    @test d.basis.monomials == [x^0, x, x^2]
    @test d.dependence == [
        MM.INDEPENDENT,
        MM.INDEPENDENT,
        MM.TRIVIAL,
    ]
    @test d.position == [MM.STANDARD, MM.STANDARD, MM.STANDARD]
end

function test_recipe(x, y, z)
    a = [x^0 * y^0]
    A = ([0], [0])
    # Corners
    c = [x, y^2]
    C = ([1, 0], [0, 2])
    ac, Ia, Ic = MB.merge_bases(b(a), b(c))
    # Dependent border
    d = [x * y * z^0]
    D = ([1], [1], [0])
    # Independent border
    e = [x * y^0 * z]
    E = ([1], [0], [1])
    f = [x^0 * y^0 * z]
    F = ([0], [0], [1])
    de, Id, Ie = MB.merge_bases(b(d), b(e))
    ae, _, _ = MB.merge_bases(b(a .* z^0), b(e))
    fe, _, _ = MB.merge_bases(b(f), b(e))
    cd, _, _ = MB.merge_bases(b(c .* z^0), b(d))
    fecd, _, Idep = MB.merge_bases(fe, cd)
    _test_recipe(
        MM.AnyDependence(FixedDependence(findall(!iszero, Ic)), ac),
        [0:1, 0:2],
        [A, C],
        ["Independent", "Dependent"],
        [:circle, :rect],
    )
    _test_recipe(
        MM.AnyDependence(FixedDependence(findall(!iszero, Ie)), de),
        [1:1, 0:1, 0:1],
        [D, E],
        ["Independent", "Dependent"],
        [:circle, :rect],
    )
    dep = MM.StaircaseDependence(FixedDependence(findall(!iszero, Idep)), fecd)
    @test sprint(show, dep) == """
StaircaseDependence for bases:
 Standard:
 MonomialBasis([z])
 Trivial Standard:
 MonomialBasis([1, y, z^2, y*z, z^3, y*z^2])
 Corners:
 MonomialBasis([x, y^2])
 Independent Border:
 MonomialBasis([x*z])
 Trivial Independent Border:
 MonomialBasis([y^2*z, x*z^2, x*y*z])
 Dependent Border:
 MonomialBasis([x*y])
"""
    @test MM.corners_basis(dep).monomials == c
    _test_recipe(
        dep,
        [0:1, 0:2, 0:3],
        [
            F,
            (zeros(Int, 6), [0, 1, 0, 1, 0, 1], [0, 0, 2, 1, 3, 2]),
            (C..., [0, 0]),
            E,
            ([0, 1, 1], [2, 0, 1], [1, 2, 1]),
            D,
        ],
        [
            "Standard",
            "Trivial Standard",
            "Corners",
            "Independent Border",
            "Trivial Independent Border",
            "Dependent Border",
        ],
        [:circle, :circle, :diamond, :rect, :rect, :rect],
    )
    return dep
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
