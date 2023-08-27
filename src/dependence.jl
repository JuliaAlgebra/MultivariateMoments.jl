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

struct StaircaseDependence
    linear::LinearDependence
    standard_or_corner::Bool
end

is_dependent(d::StaircaseDependence) = is_dependent(d.linear)
is_trivial(d::StaircaseDependence) = is_trivial(d.linear)

_is_trivial(::Bool, ::Nothing) = true
_is_trivial(a::Bool, b::Bool) = a == b

function is_standard(d::StaircaseDependence; trivial = nothing)
    return d.standard_or_corner &&
           !is_dependent(d) &&
           _is_trivial(d.linear == TRIVIAL, trivial)
end

function is_corner(d::StaircaseDependence)
    return d.standard_or_corner && is_dependent(d)
end

function _shape(d::StaircaseDependence)
    if is_standard(d)
        return :circle
    elseif is_corner(d)
        return :diamond
    else
        return :rect
    end
end

function _label(d::StaircaseDependence)
    if is_corner(d)
        return _upper("corners")
    else
        return _join(
            _label(d.linear, dependent = !is_standard(d)),
            _upper(is_standard(d) ? "standard" : "border"),
        )
    end
end

function _markercolor(d::StaircaseDependence)
    return _markercolor(d.linear)
end

"""
    struct BasisDependence{D,B<:MB.AbstractPolynomialBasis}
        basis::B
        dependence::Vector{D}
    end

The dependence between the elements of a basis.

!!! tip
    If the number of variables is 2 or 3, it can be visualized with Plots.jl.
"""
struct BasisDependence{D,B<:MB.AbstractPolynomialBasis}
    basis::B
    dependence::Vector{D}
end

function Base.isempty(d::BasisDependence)
    return isempty(d.dependence)
end

function Base.show(io::IO, d::BasisDependence)
    print(io, "BasisDependence for ")
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

function sub_basis(d::BasisDependence, I::AbstractVector{Int})
    @assert issorted(I)
    return typeof(d.basis)(d.basis.monomials[I])
end

function independent_basis(d::BasisDependence)
    return sub_basis(d, findall(!is_dependent, d.dependence))
end

function dependent_basis(d::BasisDependence)
    return sub_basis(d, findall(is_dependent, d.dependence))
end

function string_dependence(d::BasisDependence, i)
    return string(d.dependence[i])
end

function standard_basis(d::BasisDependence; kws...)
    I = findall(d.dependence) do d
        return is_standard(d; kws...)
    end
    return sub_basis(d, I)
end

function corners_basis(d::BasisDependence)
    return sub_basis(d, findall(is_corner, d.dependence))
end

function Base.convert(
    ::Type{BasisDependence{LinearDependence}},
    d::BasisDependence{StaircaseDependence},
)
    return BasisDependence(d.basis, map(d -> d.linear, d.dependence))
end

function BasisDependence{LinearDependence}(
    r,
    basis::MB.MonomialBasis{M},
) where {M}
    return BasisDependence(
        basis,
        LinearDependence[
            is_dependent!(r, i) ? DEPENDENT : INDEPENDENT for
            i in eachindex(basis.monomials)
        ],
    )
end

function is_dependent!(d::BasisDependence, i)
    return is_dependent(d.dependence[i])
end

column_compression!(::BasisDependence, ::Any) = nothing

function Base.convert(
    ::Type{BasisDependence{StaircaseDependence}},
    d::BasisDependence{LinearDependence},
)
    return BasisDependence{StaircaseDependence}(d, d.basis)
end

"""
    function BasisDependence{StaircaseDependence}(
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
function BasisDependence{StaircaseDependence}(
    r,
    basis::MB.MonomialBasis{M},
) where {M}
    if isempty(basis.monomials)
        return BasisDependence(basis, StaircaseDependence[])
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
    d = StaircaseDependence[]
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
            push!(d, StaircaseDependence(dependence(mono), true))
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
            i in eachindex(d) if !is_dependent(d[i]) for shift in vars
        ],
    )
    full_basis, I1, I2 = MB.merge_bases(full_basis, new_basis)
    deps = Vector{StaircaseDependence}(undef, length(full_basis.monomials))
    for (i, mono) in enumerate(full_basis.monomials)
        if iszero(I1[i])
            @assert !iszero(I2[i])
            if is_corner_multiple(mono, 1:(i-1), view(deps, 1:(i-1)))
                std = false
            else
                # If it was not seen before, it means it is outside the basis
                # so it is trivial standard
                @assert !is_dependent(deps[i])
                @assert is_trivial(deps[i])
                std = true
            end
            deps[i] = StaircaseDependence(dependence(mono), std)
        else
            deps[i] = d[I1[i]]
        end
    end
    return BasisDependence(full_basis, deps)
end

function _exponents(monos, i)
    return [MP.exponents(mono)[i] for mono in monos]
end

function _split_exponents(monos)
    N = MP.nvariables(monos)
    return ntuple(Base.Fix1(_exponents, monos), Val(N))
end

MP.variables(d::BasisDependence) = MP.variables(d.basis.monomials)

function MP.mindegree(d::BasisDependence, args...)
    return MP.mindegree(d.basis.monomials, args...)
end

function MP.maxdegree(d::BasisDependence, args...)
    return MP.maxdegree(d.basis.monomials, args...)
end

function _ticks(d::BasisDependence, v)
    return MP.mindegree(d, v):MP.maxdegree(d, v)
end

function _categories(d::BasisDependence{D}) where {D}
    M = eltype(d.basis.monomials)
    categories = Dict{D,Vector{M}}()
    for (i, mono) in enumerate(d.basis.monomials)
        cat = d.dependence[i]
        if !haskey(categories, cat)
            categories[cat] = M[]
        end
        push!(categories[cat], mono)
    end
    return sort!(collect(categories))
end

RecipesBase.@recipe function f(d::BasisDependence)
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
