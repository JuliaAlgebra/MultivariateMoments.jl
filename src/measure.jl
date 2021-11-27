export measure, dirac
export variables, monomials, moments

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
struct Measure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractMeasure{T}
    a::Vector{T}
    x::MVT

    function Measure{T, MT, MVT}(a::Vector{T}, x::MVT) where {T, MT, MVT}
        @assert length(a) == length(x)
        new(a, x)
    end
end
Measure(a::AbstractVector{T}, x::AbstractVector{TT}) where {T, TT <: AbstractTermLike} = Measure{T, monomialtype(TT), monovectype(x)}(monovec(a, x)...)

"""
    measure(a, X::AbstractVector{<:AbstractMonomial})

Creates a measure with moments `moment(a[i], X[i])` for each `i`.
"""
measure(a, X) = Measure(a, X)

"""
    variables(μ::AbstractMeasureLike)

Returns the variables of `μ` in decreasing order. Just like in MultivariatePolynomials, it could contain variables of zero degree in every monomial.
"""
MP.variables(μ::Measure) = variables(μ.x)

"""
    monomials(μ::AbstractMeasureLike)

Returns an iterator over the monomials of `μ` sorted in the decreasing order.
"""
MP.monomials(μ::Measure) = μ.x

"""
    moments(μ::AbstractMeasureLike)

Returns an iterator over the moments of `μ` sorted in decreasing order of monomial.
"""
moments(μ::Measure) = map((α, x) -> moment(α, x), μ.a, μ.x)

Base.:(*)(α, μ::Measure) = measure(α * μ.a, μ.x)
Base.:(*)(μ::Measure, α) = measure(μ.a * α, μ.x)
Base.:(-)(μ::Measure) = measure(-μ.a, μ.x)
function Base.:(+)(μ::Measure, ν::Measure)
    @assert μ.x == ν.x
    measure(μ.a + ν.a, μ.x)
end

"""
    dirac(X::AbstractVector{<:AbstractMoment}, s::AbstractSubstitution...)

Creates the dirac measure by evaluating the moments of `X` using `s`.

## Examples

Calling `dirac([x*y, x*y^2], x=>3, y=>2)` should the measure with moment `x*y` of value `6` and moment `x*y^2` of value `12`.
"""
function dirac(x::AbstractVector{MT}, s::MP.AbstractSubstitution...) where {MT <: AbstractMonomial}
    Measure([m(s...) for m in x], x)
end

#function truncate(μ::Measure, deg::Integer)
#    I = MP.degree.(μ.x) .<= deg
#    return Measure(μ.a[I], μ.x[I])
#end
