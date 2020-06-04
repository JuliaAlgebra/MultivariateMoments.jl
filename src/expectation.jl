const APL = AbstractPolynomialLike

function _expectation(μ::Measure{T, BT}, p::APL, f) where {T, BT}
    i = 1
    s = 0
    for (c, m) in zip(coefficients(p, BT), MB.basis_covering_monomials(BT, monomials(p)))
        while i <= length(μ.x) && m != μ.x[i]
            i += 1
        end
        if i > length(μ.x)
            error("The polynomial $p has a nonzero term $(c*m) with basis function $(m) for which the expectation is not known in $μ")
        end
        s += f(μ.a[i], c)
        i += 1
    end
    s
end

"""
    MultivariateMoments.expectation(μ::AbstractMeasureLike, p::AbstractPolynomialLike)
    MultivariateMoments.expectation(p::AbstractPolynomialLike, μ::AbstractMeasureLike)

Computes the expectation ``\\mathbb{E}_{\\mu}[p]``.
"""
function expectation end

expectation(μ::Measure, p::APL) = _expectation(μ, p, (*))
expectation(p::APL, μ::Measure) = _expectation(μ, p, (a, b) -> b * a) # a and b may be noncommutative

"""
    dot(μ::AbstractMeasureLike, p::AbstractPolynomialLike)
    dot(p::AbstractPolynomialLike, μ::AbstractMeasureLike)

See [`expectation`](@ref)
"""
LinearAlgebra.dot(μ::AbstractMeasureLike, p::APL) = expectation(μ, p)
LinearAlgebra.dot(p::APL, μ::AbstractMeasureLike) = expectation(p, μ)
