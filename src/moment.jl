export moment, moment_value, monomial

struct Moment{T,MT<:MP.AbstractMonomial} <: AbstractMoment{T}
    α::T
    x::MT
end

"""
    moment(α, m::AbstractMonomial)

Creates the moment of the monomial `m` of value `α`.
"""
moment(α, m::MP.AbstractMonomial) = Moment(α, m)

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

for f in [:variables, :nvariables, :exponents, :degree, :powers]
    @eval begin
        MP.$f(m::AbstractMomentLike) = MP.$f(MP.monomial(m))
    end
end
