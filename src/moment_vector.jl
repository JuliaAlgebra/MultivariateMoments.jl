function _check_length(values, basis)
    if length(values) != length(basis)
        throw(
            DimensionMismatch(
                "dimension must match: `values` has length `$(length(values))` and `basis` has length `$(length(basis))`",
            ),
        )
    end
end

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
struct MomentVector{T,B<:SA.AbstractBasis,V} <: AbstractMeasure{T}
    values::V
    basis::B

    function MomentVector{T,B,V}(values::V, basis::B) where {T,B,V}
        if basis isa SA.ExplicitBasis
            _check_length(values, basis)
        end
        return new{T,B,V}(values, basis)
    end
end

function moment_vector(
    values::AbstractVector{T},
    basis::SA.AbstractBasis,
) where {T}
    return MomentVector{T,typeof(basis),typeof(values)}(values, basis)
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
    return moment_vector(sorted_values, MB.SubBasis{MB.Monomial}(sorted_monos))
end

SA.basis(μ::MomentVector) = μ.basis

"""
    variables(μ::MomentVector)

Returns the variables of `μ` in decreasing order. Just like in MultivariatePolynomials, it could contain variables of zero degree in every monomial.
"""
MP.variables(μ::MomentVector) = MP.variables(SA.basis(μ))

"""
    maxdegree(μ::AbstractMeasureLike)

Returns the maximal degree of the monomials of `μ`.
"""
MP.maxdegree(μ::MomentVector) = MP.maxdegree(SA.basis(μ))

"""
    mindegree(μ::MomentVector)

Returns the minimal degree of the monomials of `μ`.
"""
MP.mindegree(μ::MomentVector) = MP.mindegree(SA.basis(μ))

"""
    extdegree(μ::MomentVector)

Returns the extremal degrees of the monomials of `μ`.
"""
MP.extdegree(μ::MomentVector) = MP.extdegree(SA.basis(μ))

"""
    moments(μ::MomentVector)

Returns an iterator over the moments of `μ` sorted in decreasing order of monomial.
"""
moments(μ::MomentVector) = map((α, x) -> moment(α, x), μ.values, SA.basis(μ))

Base.:(*)(α, μ::MomentVector) = moment_vector(α * μ.values, SA.basis(μ))
Base.:(*)(μ::MomentVector, α) = moment_vector(μ.values * α, SA.basis(μ))
Base.:(-)(μ::MomentVector) = moment_vector(-μ.values, SA.basis(μ))
function Base.:(+)(μ::MomentVector, ν::MomentVector)
    @assert SA.basis(μ) == SA.basis(ν)
    return moment_vector(μ.values + ν.values, SA.basis(μ))
end

function moment_value(μ::MomentVector, t::MP.AbstractTerm)
    return MP.coefficient(t) * moment_value(μ, MP.monomial(t))
end
function moment_value(μ::MomentVector, mono::MP.AbstractMonomial)
    return moment_value(μ, MB.Polynomial{MB.Monomial}(mono))
end

function moment_value(
    μ::MomentVector{T,<:MB.SubBasis{B}},
    p::MB.Polynomial{B},
) where {T,B}
    i = MB.monomial_index(SA.basis(μ), p.monomial)
    if isnothing(i)
        throw(ArgumentError("`$μ` does not have the moment `$p`"))
    end
    return μ.values[i]
end

function moment_value(μ::MomentVector, p::SA.AlgebraElement)
    return sum(
        coef * moment_value(μ, SA.basis(parent(p))[mono]) for
        (mono, coef) in SA.nonzero_pairs(p.coeffs)
    )
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
    return moment_vector([m(s...) for m in x], x)
end
