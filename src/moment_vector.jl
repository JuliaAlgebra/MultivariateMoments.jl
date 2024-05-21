function _check_length(values, basis)
    if SA.basis(parent) <: SA.ExplicitBasis && length(values) != length(basis)
        throw(DimensionMismatch("dimension must match: `values` has length `$(length(values))` and `basis` has length `$(length(basis))`"))
    end
end

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
struct MomentVector{T,A,V} <: AbstractMomentArray{T,A}
    values::V
    parent::A

    function MomentVector{T,A,V}(values::V, parent::A) where {T,A,V}
        if SA.basis(parent) <: SA.ExplicitBasis
            _check_length(values, SA.basis(parent))
        end
        return new{T,A,V}(values, parent)
    end
end

function moment_vector(values::AbstractVector{T}, basis::MB.SubBasis) where {T}
    parent = MB.algebra(basis)
    return MomentVector{T,typeof(parent),typeof(values)}(values, parent)
end

"""
    moment_vector(values::AbstractVector{T}, monos::AbstractVector{<:AbstractMonomial}; rtol=Base.rtoldefault(T), atol=zero(T))

Creates a measure with moments `moment(values[i], monos[i])` for each `i`.
An error is thrown if there exists `i` and `j` such that `monos[i] == monos[j]` but
`!isapprox(values[i], values[j]; rtol=rtol, atol=atol)`.
"""
function moment_vector(
    values::AbstractVector,
    monos::AbstractVector{TT};
    kws...,
) where {TT<:MP.AbstractTermLike}
    # cannot use `monomial_vector(a, x)` as it would sum the entries
    # corresponding to the same monomial.
    _check_length(values, monos)
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    sorted_values = values[σ]
    if length(monos) > length(sorted_monos)
        rev = Dict(sorted_monos[j] => j for j in eachindex(σ))
        for i in eachindex(monos)
            j = rev[monos[i]]
            if i != σ[j]
                if !isapprox(sorted_values[j], values[i]; kws...)
                    error(
                        "The monomial `$(monos[i])` occurs twice with different values: `$(values[i])` and `$(sorted_values[j])`",
                    )
                end
            end
        end
    end
    return moment_vector(values, MB.SubBasis{MB.Monomial}(sorted_monos))
end

"""
    variables(μ::AbstractMeasureLike)

Returns the variables of `μ` in decreasing order. Just like in MultivariatePolynomials, it could contain variables of zero degree in every monomial.
"""
MP.variables(μ::MomentVector) = MP.variables(μ.x)

"""
    monomials(μ::AbstractMeasureLike)

Returns an iterator over the monomials of `μ` sorted in the decreasing order.
"""
MP.monomials(μ::MomentVector) = μ.x

"""
    maxdegree(μ::AbstractMeasureLike)

Returns the maximal degree of the monomials of `μ`.
"""
MP.maxdegree(μ::MomentVector) = MP.maxdegree(MP.monomials(μ))

"""
    mindegree(μ::AbstractMeasureLike)

Returns the minimal degree of the monomials of `μ`.
"""
MP.mindegree(μ::MomentVector) = MP.mindegree(MP.monomials(μ))

"""
    extdegree(μ::AbstractMeasureLike)

Returns the extremal degrees of the monomials of `μ`.
"""
MP.extdegree(μ::MomentVector) = MP.extdegree(MP.monomials(μ))

"""
    moments(μ::AbstractMeasureLike)

Returns an iterator over the moments of `μ` sorted in decreasing order of monomial.
"""
moments(μ::MomentVector) = map((α, x) -> moment(α, x), μ.a, μ.x)

Base.:(*)(α, μ::MomentVector) = measure(α * μ.a, μ.x)
Base.:(*)(μ::MomentVector, α) = measure(μ.a * α, μ.x)
Base.:(-)(μ::MomentVector) = measure(-μ.a, μ.x)
function Base.:(+)(μ::MomentVector, ν::MomentVector)
    @assert μ.x == ν.x
    return measure(μ.a + ν.a, μ.x)
end
function _index(basis::MB.SubBasis{B}, mono) where {B}
    return get(basis, MB.Polynomial{B}(mono), nothing)
end

function moment_value(μ, mono)
    i = _index(SA.basis(μ.parent), mono)
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
