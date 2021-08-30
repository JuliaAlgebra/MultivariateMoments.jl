export measure, dirac
export variables, base_functions, moments

# If a monomial is not in basis, it does not mean that the moment is zero, it means that it is unknown/undefined
struct Measure{T, B<:MB.AbstractPolynomialBasis} <: AbstractMeasure{T}
    values::Vector{T}
    basis::B
    function Measure(values::Vector{T}, basis::B) where {T, B<:MB.AbstractPolynomialBasis}
        @assert length(values) == length(basis)
        new{T, B}(values, basis)
    end
end

basis_functions_type(::Measure{T, BT}) where {T, BT} = BT

function Measure(values::AbstractVector{T}, basis::AbstractVector{TT}) where {T, TT <: AbstractTermLike}
    b, y = monovec(values, basis)
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
MP.variables(μ::Measure) = variables(μ.basis)

"""
    monomials(μ::AbstractMeasureLike)

Returns an iterator over the monomials of `μ` sorted in the decreasing order.
"""
MP.monomials(μ::Measure) = μ.basis

"""
    base_functions(μ::AbstractMeasureLike)

Returns an iterator over the base_functions of `μ` sorted in the decreasing order.
"""
base_functions(μ::Measure) = μ.basis


"""
    moments(μ::AbstractMeasureLike)

Returns an iterator over the moments of `μ` sorted in decreasing order of monomial.
"""
moments(μ::Measure) = map((α, x) -> moment(α, x), μ.values, μ.basis)

Base.:(*)(α, μ::Measure) = measure(α * μ.values, μ.basis)
Base.:(*)(μ::Measure, α) = measure(μ.values * α, μ.basis)
Base.:(-)(μ::Measure) = measure(-μ.values, μ.basis)
function Base.:(+)(μ::Measure, ν::Measure)
    @assert μ.basis == ν.basis
    measure(μ.values + ν.values, μ.basis)
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
