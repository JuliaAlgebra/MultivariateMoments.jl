module TestQuotient

using Test
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MultivariateMoments as MM
import MutableArithmetics as MA
import SemialgebraicSets
import SemialgebraicSets as SS
import SparseArrays

b(monos) = MB.SubBasis{MB.Monomial}(MP.monomial_vector(monos))

function _test_quotient(basis, domain, expected)
    q = MM.quotient_basis(basis, domain)
    @test q isa
          MA.promote_operation(MM.quotient_basis, typeof(basis), typeof(domain))
    T = float(MP.coefficient_type(typeof(domain)))
    @test q.matrix isa SparseArrays.SparseMatrixCSC{T,Int}
    @test q.basis == basis
    # The Macaulay matrix path is numerical so the entries may be rounded
    @test q.matrix ≈ expected
    return
end

function test_groebner(x, y)
    # A single equation is a Gröbner basis
    domain = SS.@set x^2 + 1 == 0
    @test SS.is_gröbner_basis(SS.ideal(domain))
    _test_quotient(b([1, x, x^2, x^3]), domain, [1 0 -1 0; 0 1 0 -1])
    _test_quotient(b(typeof(x * y)[]), domain, zeros(0, 0))
    return
end

function test_fixed_variables(x, y)
    domain = SS.@set x == 1
    _test_quotient(b([1, y, x, x * y]), domain, [1 0 1 0; 0 1 0 1])
    _test_quotient(b([1, y, x, x * y]), SS.ideal(domain), [1 0 1 0; 0 1 0 1])
    return
end

function test_macaulay(x, y)
    # The equation is redundant so that it is not a Gröbner basis
    domain = SS.@set x^2 + 1.0 == 0 && 2x^2 + 2.0 == 0
    @test !SS.is_gröbner_basis(SS.ideal(domain))
    # `y` is not a column of the Macaulay matrix of degree 2
    _test_quotient(b([1, y, x^2]), domain, [1 0 -1; 0 1 0])
    # Same matrix as the remainders modulo the Gröbner basis
    _test_quotient(b([1, x, x^2, x^3]), domain, [1 0 -1 0; 0 1 0 -1])
    # The Macaulay matrix of degree 1 has no rows
    _test_quotient(b([1, x]), domain, [1 0; 0 1])
    _test_quotient(b(typeof(x * y)[]), domain, zeros(0, 0))
    @test !SS.is_gröbner_basis(SS.ideal(domain))
    return
end

function test_moment_vector(x, y)
    domain = SS.@set x^2 + 1 == 0
    q = MM.quotient_basis(b([1, x, x^2, x^3]), domain)
    ν = MM.moment_vector([2.0, 3.0], q)
    @test ν.basis == q.basis
    @test ν.values == [2, 3, -2, -3]
    @test_throws DimensionMismatch MM.moment_vector([1.0], q)
    return
end

function _in_kernel(q, p)
    a = MP.coefficients(p, MB.keys_as_monomials(q.basis))
    return isapprox(q.matrix * a, zeros(size(q.matrix, 1)), atol = 1e-12)
end

function test_degree_fall(x, y)
    # `y * (x * y - 1) - x * (y^2 - 1) = x - y` is in the ideal but it is not
    # in its truncation at degree 1 or 2
    domain = SS.@set x * y == 1.0 && y^2 == 1.0
    linear = b(MP.monomials((x, y), 0:1))
    # By default, the ideal is truncated at the maximal degree of the basis
    for d in 1:2
        @test !_in_kernel(
            MM.quotient_basis(b(MP.monomials((x, y), 0:d)), domain),
            x - y,
        )
        @test !_in_kernel(
            MM.quotient_basis(linear, SS.truncated(domain, d)),
            x - y,
        )
    end
    @test _in_kernel(
        MM.quotient_basis(b(MP.monomials((x, y), 0:3)), domain),
        x - y,
    )
    # The truncation degree can be larger than the degree of the basis
    q = MM.quotient_basis(linear, SS.truncated(domain, 3))
    @test q isa MA.promote_operation(
        MM.quotient_basis,
        typeof(linear),
        typeof(SS.truncated(domain, 3)),
    )
    @test _in_kernel(q, x - y)
    @test !_in_kernel(q, x)
    @test !SS.is_gröbner_basis(SS.ideal(domain))
    return
end

function test_truncated_smaller_degree(x, y)
    # The truncation is used even if the generator is a Gröbner basis
    domain = SS.truncated(SS.@set(x^2 + 1.0 == 0), 2)
    # `x^3` is not in the Macaulay matrix of degree 2 so it is independent
    _test_quotient(b([1, x, x^2, x^3]), domain, [1 0 -1 0; 0 1 0 0; 0 0 0 1])
    return
end

function _dependent(I, vars)
    d = MM.BasisDependence{MM.LinearDependence}(I, vars)
    return MB.keys_as_monomials(MM.dependent_basis(d))
end

function test_dependence(x, y)
    I = SS.truncated(SS.ideal(SS.@set(x^2 + 1.0 == 0)), 2)
    d = MM.BasisDependence{MM.LinearDependence}(I, MP.variables(x))
    # `x` is not a column of the Macaulay matrix
    @test MB.keys_as_monomials(d.basis) == [1, x^2]
    @test _dependent(I, MP.variables(x)) == [x^2]
    @test _dependent(I, MP.variables(x * y)) == [x^2]
    # `x - y` is in the truncation at degree 3 but not 2
    domain = SS.@set x * y == 1.0 && y^2 == 1.0
    lead = MP.leading_monomial(x - y)
    for d in 2:3
        dep = _dependent(SS.truncated(SS.ideal(domain), d), MP.variables(x * y))
        @test (lead in dep) == (d == 3)
        @test y^2 in dep
    end
    # The Macaulay matrix of degree 1 has no rows
    @test isempty(
        _dependent(SS.truncated(SS.ideal(domain), 1), MP.variables(x * y)),
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
    DynamicPolynomials.@polyvar x y
    TestQuotient.runtests(x, y)
end

import TypedPolynomials
@testset "TypedPolynomials" begin
    TypedPolynomials.@polyvar x y
    TestQuotient.runtests(x, y)
end
