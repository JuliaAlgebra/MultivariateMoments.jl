export Moment, moment, value

struct Moment{T, MT <: AbstractMonomial} <: AbstractMoment{T}
    α::T
    x::MT
end

"""
    moment(α, m::AbstractMonomial)

Creates the moment of the monomial `m` and of value `α`.
"""
moment(α, m::AbstractMonomial) = Moment(α, m)

"""
    value(m::AbstractMomentLike{T})::T where T

Returns the value of the moment `m`.

## Examples

Calling `value(moment(3.1, x*y^2))` should return `3.1`.
"""
value(m::Moment) = m.α

MP.monomial(m::Moment) = m.x

for f in [:variables, :nvariables, :exponents, :degree, :powers]
    @eval begin
        MP.$f(m::AbstractMomentLike) = $f(monomial(m))
    end
end
