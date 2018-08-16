const APL = AbstractPolynomialLike

function _expectation(μ::Measure, p::APL, f)
    i = 1
    s = 0
    for t in terms(p)
        while i <= length(μ.x) && monomial(t) != μ.x[i]
            i += 1
        end
        if i > length(μ.x)
            error("The polynomial $p has a nonzero term $t with monomial $(t.x) for which the expectation is not known in $μ")
        end
        s += f(μ.a[i], coefficient(t))
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
Compat.LinearAlgebra.dot(μ::AbstractMeasureLike, p::APL) = expectation(μ, p)
Compat.LinearAlgebra.dot(p::APL, μ::AbstractMeasureLike) = expectation(p, μ)
