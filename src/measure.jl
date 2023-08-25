# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
struct Measure{T,MT<:MP.AbstractMonomial,MVT<:AbstractVector{MT}} <:
       AbstractMeasure{T}
    a::Vector{T}
    x::MVT

    function Measure{T,MT,MVT}(a::Vector{T}, x::MVT) where {T,MT,MVT}
        @assert length(a) == length(x)
        return new(a, x)
    end
end
function Measure(
    a::AbstractVector{T},
    x::AbstractVector{TT};
    kws...,
) where {T,TT<:MP.AbstractTermLike}
    # cannot use `monomial_vector(a, x)` as it would sum the entries
    # corresponding to the same monomial.
    if length(a) != length(x)
        throw(
            ArgumentError("There should be as many coefficient than monomials"),
        )
    end
    σ, X = MP.sort_monomial_vector(x)
    b = a[σ]
    if length(x) > length(X)
        rev = Dict(X[j] => j for j in eachindex(σ))
        for i in eachindex(x)
            j = rev[x[i]]
            if i != σ[j]
                if !isapprox(b[j], a[i]; kws...)
                    @warn(
                        "The monomial `$(x[i])` occurs twice with different values: `$(a[i])` and `$(b[j])`, keeping the first one.",
                    )
                end
            end
        end
    end
    return Measure{T,MP.monomial_type(TT),typeof(X)}(b, X)
end

"""
    measure(a::AbstractVector{T}, X::AbstractVector{<:AbstractMonomial}; rtol=Base.rtoldefault(T), atol=zero(T))

Creates a measure with moments `moment(a[i], X[i])` for each `i`.
An error is thrown if there exists `i` and `j` such that `X[i] == X[j]` but
`!isapprox(a[i], a[j]; rtol=rtol, atol=atol)`.
"""
measure(a, X; kws...) = Measure(a, X; kws...)
function measure(a, basis::MB.MonomialBasis; kws...)
    return measure(a, basis.monomials; kws...)
end

"""
    variables(μ::AbstractMeasureLike)

Returns the variables of `μ` in decreasing order. Just like in MultivariatePolynomials, it could contain variables of zero degree in every monomial.
"""
MP.variables(μ::Measure) = MP.variables(μ.x)

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
    return measure(μ.a + ν.a, μ.x)
end

# TODO Move to MultivariateBases
function _monomial_index(monos::AbstractVector, mono)
    i = searchsortedfirst(monos, mono)
    if i in eachindex(monos) && mono == monos[i]
        return i
    else
        return
    end
end

function _index(basis::MB.MonomialBasis, mono)
    return _monomial_index(basis.monomials, mono)
end

function moment_value(μ, mono)
    i = _monomial_index(μ.x, mono)
    if isnothing(i)
        throw(ArgumentError("`$μ` does not have the moment `$mono`"))
    end
    return μ.a[i]
end

"""
    dirac(X::AbstractVector{<:AbstractMoment}, s::AbstractSubstitution...)

Creates the dirac measure by evaluating the moments of `X` using `s`.

## Examples

Calling `dirac([x*y, x*y^2], x=>3, y=>2)` should the measure with moment `x*y` of value `6` and moment `x*y^2` of value `12`.
"""
function dirac(
    x::AbstractVector{MT},
    s::MP.AbstractSubstitution...,
) where {MT<:MP.AbstractMonomial}
    return Measure([m(s...) for m in x], x)
end
