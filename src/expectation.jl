const APL = AbstractPolynomialLike

function _expectation(μ::Measure{T, BT}, p::APL, f) where {T, BT}
    i = 1
    s = 0
    basis = MB.basis_covering_monomials(BT, monomials(p))
    for (c, m) in zip(coefficients(p, basis), basis)
        while i <= length(μ.basis) && m != μ.basis[i]
            i += 1
        end
        if i > length(μ.basis)
            error("The base function $m has a non-zero multiplier $c in $p, 
                  but its expectation is not known in $μ")
        end
        s += f(μ.values[i], c)
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
