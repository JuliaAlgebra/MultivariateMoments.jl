struct Moment{T,P} <: AbstractMoment{T}
    α::T
    polynomial::P
end

"""
    moment(α, p)

Creates the moment of the polynomial `p` of value `α`.
"""
moment(α, p) = Moment(α, p)

"""
    moment_value(m::AbstractMomentLike)

Returns the value of the moment `m`.

## Examples

Calling `moment_value(moment(3.1, x*y^2))` should return `3.1`.
"""
moment_value(m::Moment) = m.α

for f in [:variables, :nvariables]
    @eval begin
        MP.$f(m::AbstractMomentLike) = MP.$f(m.polynomial)
    end
end
