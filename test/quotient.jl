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
    @test q.matrix == expected
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

function test_degree_fall(x, y)
    # `y * (x * y - 1) - x * (y^2 - 1) = x - y` is in the ideal but it is not
    # in its truncation at degree 1 or 2
    domain = SS.@set x * y == 1.0 && y^2 == 1.0
    for d in 1:2
        q = MM.quotient_basis(b(MP.monomials((x, y), 0:d)), domain)
        @test !iszero(
            q.matrix * MP.coefficients(x - y, MB.keys_as_monomials(q.basis)),
        )
    end
    q = MM.quotient_basis(b(MP.monomials((x, y), 0:3)), domain)
    a = MP.coefficients(x - y, MB.keys_as_monomials(q.basis))
    @test q.matrix * a ≈ zeros(size(q.matrix, 1)) atol = 1e-12
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
