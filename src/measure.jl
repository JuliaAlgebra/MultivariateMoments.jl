export measure, dirac
export variables, base_functions, moments

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
struct Measure{T, B<:MB.AbstractPolynomialBasis} <: AbstractMeasure{T}
    a::Vector{T}
    x::B
    function Measure(a::Vector{T}, x::B) where {T, B<:MB.AbstractPolynomialBasis}
        @assert length(a) == length(x)
        new{T, B}(a, x)
    end
end

basis_functions_type(::Measure{T, BT}) where {T,BT} = BT

function Measure(a::AbstractVector{T}, x::AbstractVector{TT}) where {T, TT <: AbstractTermLike}
    b, y = monovec(a, x)
    return Measure(b, MB.MonomialBasis(y))
end

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
    base_functions(μ::AbstractMeasureLike)

Returns an iterator over the base_functions of `μ` sorted in the decreasing order.
"""
base_functions(μ::Measure) = μ.x


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
