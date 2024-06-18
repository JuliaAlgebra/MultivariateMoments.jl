function _expectation(
    μ::MomentVector{S,<:MB.SubBasis{B}},
    p::SA.AlgebraElement{<:MB.Algebra{BT,B},T},
    f,
) where {S,T,BT,B}
    i = firstindex(μ.basis)
    s = zero(MA.promote_operation(*, S, T))
    for (k, v) in SA.nonzero_pairs(SA.coeffs(p))
        mono = SA.basis(p)[k]
        while i in eachindex(μ.basis) && mono != μ.basis[i]
            i += 1
        end
        if !(i in eachindex(μ.basis))
            error(
                "The polynomial $p has a nonzero term $mono with coefficient $v for which the expectation is not known in $μ",
            )
        end
        s += f(μ.values[i], v)
        i += 1
    end
    return s
end

function _expectation(
    μ::MomentVector{S,<:MB.SubBasis{B}},
    p::SA.AlgebraElement{<:MB.Algebra},
    f,
) where {S,B}
    basis = MB.FullBasis{B,MP.monomial_type(typeof(p))}()
    return _expectation(μ, MB.algebra_element(SA.coeffs(p, basis), basis), f)
end

function _expectation(μ::MomentVector, p::MP.AbstractPolynomialLike, f)
    return _expectation(
        μ,
        MB.algebra_element(
            MB.sparse_coefficients(MP.polynomial(p)),
            MB.FullBasis{MB.Monomial,MP.monomial_type(p)}(),
        ),
        f,
    )
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
LinearAlgebra.dot(μ::AbstractMeasureLike, p::MP.AbstractPolynomialLike) =
    expectation(μ, p)
function LinearAlgebra.dot(p::MP.AbstractPolynomialLike, μ::AbstractMeasureLike)
    return expectation(p, μ)
end
# Need to split it and not use `_APL` to avoid ambiguity
function LinearAlgebra.dot(μ::AbstractMeasureLike, p::SA.AlgebraElement)
    return expectation(μ, p)
end
function LinearAlgebra.dot(p::SA.AlgebraElement, μ::AbstractMeasureLike)
    return expectation(p, μ)
end
