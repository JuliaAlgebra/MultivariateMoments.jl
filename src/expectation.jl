function _expectation(μ::MomentVector{S,<:MB.SubBasis{MB.Monomial}}, p::_APL{T}, f) where {S,T}
    i = firstindex(μ.basis)
    s = zero(MA.promote_operation(*, S, T))
    for t in MP.terms(p)
        while i in eachindex(μ.basis) && MP.monomial(t) != μ.basis.monomials[i]
            i += 1
        end
        if !(i in eachindex(μ.basis))
            error(
                "The polynomial $p has a nonzero term $t with monomial $(MP.monomial(t)) for which the expectation is not known in $μ",
            )
        end
        s += f(μ.values[i], MP.coefficient(t))
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

expectation(μ::MomentVector, p::_APL) = _expectation(μ, p, (*))
expectation(p::_APL, μ::MomentVector) = _expectation(μ, p, (a, b) -> b * a) # a and b may be noncommutative

"""
    dot(μ::AbstractMeasureLike, p::AbstractPolynomialLike)
    dot(p::AbstractPolynomialLike, μ::AbstractMeasureLike)

See [`expectation`](@ref)
"""
LinearAlgebra.dot(μ::AbstractMeasureLike, p::_APL) = expectation(μ, p)
LinearAlgebra.dot(p::_APL, μ::AbstractMeasureLike) = expectation(p, μ)
