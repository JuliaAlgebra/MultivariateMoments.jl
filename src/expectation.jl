const _AE{B,T} = SA.AlgebraElement{T,<:SA.StarAlgebra{<:MB.Variables{B}}}

function SA.promote_with_map(μ::MomentVector, basis, map)
    return moment_vector(μ.values, basis), map
end

function SA.promote_basis_with_maps(
    μ::MomentVector,
    p::Union{SA.AlgebraElement,SA.AbstractBasis,SA.AbstractStarAlgebra},
)
    _p, _μ = SA.promote_basis_with_maps(p, SA.basis(μ))
    return SA.maybe_promote(μ, _μ...), _p
end

function _same_basis_expectation(
    μ::MomentVector{S,<:MB.SubBasis{B}},
    p::_AE{B},
    f,
) where {S,B}
    s = zero(MA.promote_operation(*, S, eltype(p)))
    for (k, v) in SA.nonzero_pairs(SA.coeffs(p))
        el = SA.basis(p)[k]
        i = get(μ.basis, el, nothing)
        if isnothing(i)
            error(
                "The polynomial $p has a nonzero term $el with coefficient $v for which the expectation is not known in $μ",
            )
        end
        s += f(μ.values[i], v)
    end
    return s
end

function _expectation(
    μ::MomentVector{S,<:MB.SubBasis{B}},
    p::_AE{B},
    f,
) where {S,B}
    return _same_basis_expectation(SA.promote_basis(μ, p)..., f)
end

function _expectation(
    μ::MomentVector{S,<:MB.SubBasis{B}},
    p::_AE,
    f,
) where {S,B}
    basis = MB.FullBasis{B}(MP.variables(p))
    return _expectation(μ, MB.algebra_element(SA.coeffs(p, basis), basis), f)
end

function _expectation(μ::MomentVector, p::MP.AbstractPolynomialLike, f)
    return _expectation(μ, MB.algebra_element(p), f)
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
function LinearAlgebra.dot(μ::AbstractMeasureLike, p::MP.AbstractPolynomialLike)
    return expectation(μ, p)
end
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
