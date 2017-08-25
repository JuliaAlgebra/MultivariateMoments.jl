export Measure, dirac

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
struct Measure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractMeasure{T}
    a::Vector{T}
    x::MVT

    function Measure{T, MT, MVT}(a::Vector{T}, x::MVT) where {T, MT, MVT}
        @assert length(a) == length(x)
        new(a, x)
    end
end

measure(a, X) = Measure(a, X)
(*)(α, μ::Measure) = measure(α * μ.a, μ.x)
(*)(μ::Measure, α) = measure(μ.a * α, μ.x)

function (+)(μ::Measure, ν::Measure)
    @assert μ.x == ν.x
    measure(μ.a + ν.a, μ.x)
end

MP.monomials(μ::Measure) = μ.x
moments(μ::Measure) = map((α, x) -> moment(α, x), μ.a, μ.x)

Measure(a::Vector{T}, x::AbstractVector{TT}) where {T, TT <: AbstractTermLike} = Measure{T, monomialtype(TT), monovectype(x)}(monovec(a, x)...)

"""
    dirac(X::AbstractVector{<:AbstractMoment}, s::AbstractSubstitution...)

Creates the dirac measure by evaluating the moments of `X` using `s`.

## Examples

Calling `dirac([x*y, x*y^2], x=>3, y=>2)` should the measure with moment `x*y` of value `6` and moment `x*y^2` of value `12`.
"""
function dirac(x::AbstractVector{MT}, s::MP.AbstractSubstitution...) where {MT <: AbstractMonomial}
    Measure([m(s...) for m in x], x)
end
