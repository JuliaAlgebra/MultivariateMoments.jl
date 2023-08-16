import MultivariatePolynomials as MP
import MultivariateBases as MB
import RecipesBase
import Colors

"""
    @enum LinearDependence INDEPENDENT TRIVIAL DEPENDENT

Linear dependence of the element of a basis representing the indices of the rows
of a [`MacaulayNullspace`]. `DEPENDENT` indicates that it is linearly dependent
to rows that appear earlier in the basis. `TRIVIAL` indicates that the element
was not in the original basis so it is trivially independent.
"""
@enum LinearDependence INDEPENDENT TRIVIAL DEPENDENT

is_dependent(d::LinearDependence) = d == DEPENDENT
is_trivial(d::LinearDependence) = d == TRIVIAL

_shape(d::LinearDependence) = is_dependent(d) ? :rect : :circle

function _upper(a)
    if isempty(a)
        return a
    else
        return uppercase(a[1]) * a[2:end]
    end
end

function __join(a, b)
    if isempty(a) || isempty(b)
        return a * b
    else
        return a * " " * b
    end
end

_join(a, b) = __join(a, _upper(b))

function _label(d::LinearDependence; dependent::Bool = true)
    s = ""
    s = _join(s, is_trivial(d) ? "trivial" : "")
    if dependent
        s = _join(s, (is_dependent(d) ? "" : "in") * "dependent")
    end
    return s
end

function _markercolor(d::LinearDependence)
    if is_trivial(d)
        return Colors.JULIA_LOGO_COLORS.purple
    else
        if is_dependent(d)
            return Colors.JULIA_LOGO_COLORS.green
        else
            return Colors.JULIA_LOGO_COLORS.blue
        end
    end
end

@enum StaircasePosition STANDARD CORNER BORDER

function _shape(d::Tuple{StaircasePosition,LinearDependence})
    if d[1] == STANDARD
        return :circle
    elseif d[1] == CORNER
        return :diamond
    else
        @assert d[1] == BORDER
        return :rect
    end
end

function _label(d::Tuple{StaircasePosition,LinearDependence})
    if d[1] == CORNER
        return _upper("corners")
    else
        return _join(
            _label(d[2], dependent = d[1] != STANDARD),
            lowercase(string(d[1])),
        )
    end
end

function _markercolor(d::Tuple{StaircasePosition,LinearDependence})
    return _markercolor(d[2])
end

abstract type AbstractDependence end

function Base.isempty(d::AbstractDependence)
    return isempty(d.dependence)
end

function Base.show(io::IO, d::AbstractDependence)
    print(io, "$(nameof(typeof(d))) for ")
    if isempty(d)
        println(io, "an empty basis")
    else
        println(io, "bases:")
    end
    for (cat, monos) in _categories(d)
        println(io, " ", _label(cat), ":")
        println(io, " ", MB.MonomialBasis(monos))
    end
    return
end

function sub_basis(d::AbstractDependence, I::AbstractVector{Int})
    @assert issorted(I)
    return typeof(d.basis)(d.basis.monomials[I])
end

function independent_basis(d::AbstractDependence)
    return sub_basis(d, findall(d -> !is_dependent(d), d.dependence))
end

function dependent_basis(d::AbstractDependence)
    return sub_basis(d, findall(d -> is_dependent(d), d.dependence))
end

"""
    struct AnyDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
        basis::B
        dependence::Vector{LinearDependence}
    end

The independent and dependent can be arbitrary disjoint bases.

!!! tip
    If the number of variables is 2 or 3, it can be visualized with Plots.jl.
"""
struct AnyDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
    basis::B
    dependence::Vector{LinearDependence}
end

function string_dependence(d::AnyDependence, i)
    return string(d.dependence[i])
end

_category_type(::AnyDependence) = LinearDependence
_category(d::AnyDependence, i) = d.dependence[i]

"""
    struct StaircaseDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
        basis::B
        dependence::Vector{LinearDependence}
        position::Vector{StaircasePosition}
    end

No independent is a multiple of a dependent and each dependent can be obtained
by multiplying an independent with a variable.

!!! tip
    If the number of variables is 2 or 3, it can be visualized with Plots.jl.
"""
struct StaircaseDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
    basis::B
    dependence::Vector{LinearDependence}
    position::Vector{StaircasePosition}
end

_in_basis(::Bool, ::Nothing) = true
_in_basis(a::Bool, b::Bool) = a == b

function standard_basis(d::AbstractDependence; in_basis = nothing)
    I = findall(eachindex(d.position)) do i
        return d.position[i] == STANDARD &&
               _in_basis(!is_trivial(d.dependence[i]), in_basis)
    end
    return sub_basis(d, I)
end

function corners_basis(d::AbstractDependence)
    return sub_basis(d, findall(isequal(CORNER), d.position))
end

function string_dependence(d::StaircaseDependence, i)
    return string(d.position[i]) * " " * string(d.dependence[i])
end

function _category_type(::StaircaseDependence)
    return Tuple{StaircasePosition,LinearDependence}
end
_category(d::StaircaseDependence, i) = d.position[i], d.dependence[i]

function Base.convert(::Type{AnyDependence}, d::StaircaseDependence)
    return AnyDependence(d.basis, d.dependence)
end

function AnyDependence(r, basis::MB.MonomialBasis{M}) where {M}
    return AnyDependence(
        basis,
        LinearDependence[
            is_dependent!(r, i) ? DEPENDENT : INDEPENDENT for
            i in eachindex(basis.monomials)
        ],
    )
end

function is_dependent!(d::AbstractDependence, i)
    return is_dependent(d.dependence[i])
end

column_compression!(::AbstractDependence, ::Any) = nothing

function Base.convert(::Type{StaircaseDependence}, d::AnyDependence)
    return StaircaseDependence(d, d.basis)
end

"""
    function StaircaseDependence(
        is_dependent::Function,
        basis::MB.AbstractPolynomialBasis,
    )

Computes the set of standard monomials using the *greedy sieve* algorithm
presented in [LLR08, Algorithm 1]. A monomial outside of `basis` is assumed
to be independent. Otherwise, if its index in `basis` is `i`, `is_dependent`
returns whether it is dependent.

[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp.
*Semidefinite characterization and computation of zero-dimensional real radical ideals.*
Foundations of Computational Mathematics 8 (2008): 607-647.
"""
function StaircaseDependence(r, basis::MB.MonomialBasis{M}) where {M}
    if isempty(basis.monomials)
        return StaircaseDependence(
            basis,
            LinearDependence[],
            StaircasePosition[],
        )
    end
    function dependence(mono)
        i = _index(basis, mono)
        return if isnothing(i)
            TRIVIAL
        else
            is_dependent!(r, i) ? DEPENDENT : INDEPENDENT
        end
    end
    vars = MP.variables(basis)
    full_basis =
        MB.maxdegree_basis(typeof(basis), vars, MP.maxdegree(basis.monomials))
    d = LinearDependence[]
    s = StaircasePosition[]
    # This sieve of [LLR08, Algorithm 1] is a performance improvement but not only.
    # It also ensures that the standard monomials have the "staircase structure".
    function is_corner_multiple(mono, indices, dependence)
        for i in eachindex(dependence)
            if is_dependent(dependence[i]) &&
               MP.divides(full_basis.monomials[indices[i]], mono)
                return true
            end
        end
        return false
    end
    keep = Int[]
    # Compute standard monomials and corners
    for (i, mono) in enumerate(full_basis.monomials)
        if !is_corner_multiple(mono, keep, d)
            push!(keep, i)
            push!(d, dependence(mono))
            push!(s, is_dependent(d[end]) ? CORNER : STANDARD)
        end
    end
    column_compression!(
        r,
        Int[keep[i] for i in eachindex(d) if !is_dependent(d[i])],
    )
    full_basis = typeof(full_basis)(full_basis.monomials[keep])
    new_basis = MB.MonomialBasis(
        eltype(basis.monomials)[
            full_basis.monomials[i] * shift for
            i in eachindex(s) if !is_dependent(d[i]) for shift in vars
        ],
    )
    full_basis, I1, I2 = MB.merge_bases(full_basis, new_basis)
    deps = Vector{LinearDependence}(undef, length(full_basis.monomials))
    position = Vector{StaircasePosition}(undef, length(full_basis.monomials))
    for (i, mono) in enumerate(full_basis.monomials)
        if iszero(I1[i])
            deps[i] = dependence(mono)
            @assert !iszero(I2[i])
            if is_corner_multiple(mono, 1:(i-1), view(deps, 1:(i-1)))
                position[i] = BORDER
            else
                # If it was not seen before, it means it is outside the basis
                # so it is trivial standard
                @assert !is_dependent(deps[i])
                @assert is_trivial(deps[i])
                position[i] = STANDARD
            end
        else
            deps[i] = d[I1[i]]
            position[i] = s[I1[i]]
        end
    end
    return StaircaseDependence(full_basis, deps, position)
end

function _exponents(monos, i)
    return [MP.exponents(mono)[i] for mono in monos]
end

function _split_exponents(monos)
    N = MP.nvariables(monos)
    return ntuple(Base.Fix1(_exponents, monos), Val(N))
end

MP.variables(d::AbstractDependence) = MP.variables(d.basis.monomials)

function MP.mindegree(d::AbstractDependence, args...)
    return MP.mindegree(d.basis.monomials, args...)
end

function MP.maxdegree(d::AbstractDependence, args...)
    return MP.maxdegree(d.basis.monomials, args...)
end

function _ticks(d::AbstractDependence, v)
    return MP.mindegree(d, v):MP.maxdegree(d, v)
end

function _categories(d::AbstractDependence)
    M = eltype(d.basis.monomials)
    categories = Dict{_category_type(d),Vector{M}}()
    for (i, mono) in enumerate(d.basis.monomials)
        cat = _category(d, i)
        if !haskey(categories, cat)
            categories[cat] = M[]
        end
        push!(categories[cat], mono)
    end
    return sort!(collect(categories))
end

RecipesBase.@recipe function f(d::AbstractDependence)
    vars = MP.variables(d)
    t = _ticks.(Ref(d), vars)
    aspect_ratio --> :equal # defaults to `:auto`
    xticks --> t[1]
    yticks --> t[2]
    if length(t) >= 3
        zticks --> t[3]
    end
    for (cat, monos) in _categories(d)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> _shape(cat)
            label := _label(cat)
            _split_exponents(monos)
        end
    end
end
