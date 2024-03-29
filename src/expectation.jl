function _expectation(μ::Measure{S}, p::_APL{T}, f) where {S,T}
    i = 1
    s = zero(MA.promote_operation(*, S, T))
    for t in MP.terms(p)
        while i <= length(μ.x) && MP.monomial(t) != μ.x[i]
            i += 1
        end
        if i > length(μ.x)
            error(
                "The polynomial $p has a nonzero term $t with monomial $(t.x) for which the expectation is not known in $μ",
            )
        end
        s += f(μ.a[i], MP.coefficient(t))
        i += 1
    end
    return s
end

"""
    MultivariateMoments.expectation(μ::AbstractMeasureLike, p::AbstractPolynomialLike)
    MultivariateMoments.expectation(p::AbstractPolynomialLike, μ::AbstractMeasureLike)

Computes the expectation ``\\mathbb{E}_{\\mu}[p]``.
"""
function expectation end

expectation(μ::Measure, p::_APL) = _expectation(μ, p, (*))
expectation(p::_APL, μ::Measure) = _expectation(μ, p, (a, b) -> b * a) # a and b may be noncommutative

"""
    dot(μ::AbstractMeasureLike, p::AbstractPolynomialLike)
    dot(p::AbstractPolynomialLike, μ::AbstractMeasureLike)

See [`expectation`](@ref)
"""
LinearAlgebra.dot(μ::AbstractMeasureLike, p::_APL) = expectation(μ, p)
LinearAlgebra.dot(p::_APL, μ::AbstractMeasureLike) = expectation(p, μ)
