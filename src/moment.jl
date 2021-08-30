export moment, moment_value, monomial, base_function

struct Moment{T, PT <: AbstractPolynomialLike} <: AbstractMoment{T}
    α::T
    x::PT
end

"""
    moment(α, m::AbstractPolynomialLike)

Creates the moment of the base function `m` of value `α`.
"""
moment(α, m::AbstractPolynomialLike) = Moment(α, m)

"""
    moment_value(m::AbstractMomentLike)

Returns the value of the moment `m`.

## Examples

Calling `moment_value(moment(3.1, x*y^2))` should return `3.1`.
"""
moment_value(m::Moment) = m.α

"""
    monomial(m::AbstractMomentLike)

Returns the monomial of the moment `m`.

## Examples

Calling `monomial(moment(3.1, x*y^2))` should return `x*y^2`.
"""
MP.monomial(m::Moment) = m.x

"""
    base_function(m::AbstractMomentLike)

Returns the base_function of the moment `m`.

## Examples

Calling `base_function(moment(3.1, x*y^2))` should return `x*y^2`.
"""
base_function(m::Moment) = m.x

for f in [:variables, :nvariables, :exponents, :powers]
    @eval begin
        MP.$f(m::AbstractMomentLike) = $f(base_function(m))
    end
end

MP.degree(m::AbstractMomentLike) = MP.maxdegree(base_function(m))
