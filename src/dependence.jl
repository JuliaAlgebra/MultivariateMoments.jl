import MultivariatePolynomials as MP
import MultivariateBases as MB
import RecipesBase

abstract type AbstractDependence end

"""
    struct AnyDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
        independent::B
        dependent::B
    end

The independent and dependent can be arbitrary disjoint bases.

!!! tip
    If the number of variables is 2 or 3, it can be visualized with Plots.jl.
"""
struct AnyDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
    independent::B
    dependent::B
end

function Base.isempty(d::AnyDependence)
    return isempty(d.independent.monomials) && isempty(d.dependent.monomials)
end

function Base.show(io::IO, d::AnyDependence)
    println(io, "AnyDependence")
    println(io, "with independent basis:")
    println(io, d.independent)
    println(io, "and dependent basis:")
    println(io, d.dependent)
    return
end

"""
    struct StaircaseDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
        standard::B
        corners::B
        dependent_border::B
        independent_border::B
    end

No independent is a multiple of a dependent and each dependent can be obtained
by multiplying an independent with a variable.

!!! tip
    If the number of variables is 2 or 3, it can be visualized with Plots.jl.
"""
struct StaircaseDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
    standard::B
    corners::B
    dependent_border::B
    independent_border::B
end

function Base.isempty(d::StaircaseDependence)
    return isempty(d.standard.monomials) &&
           isempty(d.corners.monomials) &&
           isempty(d.dependent_border.monomials) &&
           isempty(d.independent_border.monomials)
end

function Base.show(io::IO, d::StaircaseDependence)
    println(io, "StaircaseDependence")
    println(io, "with standard basis:")
    println(io, d.standard)
    println(io, "and corners basis:")
    println(io, d.corners)
    println(io, "and dependent border basis:")
    println(io, d.dependent_border)
    println(io, "and independent border basis:")
    print(io, d.independent_border)
    return
end

function Base.convert(::Type{AnyDependence}, d::StaircaseDependence)
    dependent = MB.merge_bases(d.corners, d.dependent_border)[1]
    return AnyDependence(d.standard, dependent)
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
function StaircaseDependence(
    is_dependent::Function,
    basis::MB.MonomialBasis{M},
) where {M}
    standard = M[]
    corners = M[]
    vars = MP.variables(basis)
    for mono in MP.monomials(vars, 0:MP.maxdegree(basis.monomials))
        # This sieve of [LLR08, Algorithm 1] is a performance improvement but not only.
        # It also ensures that the standard monomials have the "staircase structure".
        if !any(Base.Fix2(MP.divides, mono), corners)
            i = _index(basis, mono)
            if isnothing(i) || is_dependent(i)
                push!(standard, mono)
            else
                push!(corners, mono)
            end
        end
    end
    dependent_border = M[]
    independent_border = M[]
    for mono in standard
        for shift in vars
            border = shift * mono
            if isnothing(_monomial_index(standard, border)) &&
               isnothing(_monomial_index(corners, border))
                i = _index(basis, mono)
                if isnothing(i) || is_dependent(i)
                    push!(independent_border, border)
                else
                    push!(dependent_border, border)
                end
            end
        end
    end
    return StaircaseDependence(
        MB.MonomialBasis(standard),
        MB.MonomialBasis(corners),
        MB.MonomialBasis(dependent_border),
        MB.MonomialBasis(independent_border),
    )
end

function _exponents(monos, i)
    return [MP.exponents(mono)[i] for mono in monos]
end

function _split_exponents(monos)
    N = MP.nvariables(monos)
    return ntuple(Base.Fix1(_exponents, monos), Val(N))
end

function _split_exponents(basis::MB.AbstractPolynomialBasis)
    return _split_exponents(basis.monomials)
end

function _deg(deg, b, args...)
    if isempty(b.monomials)
        return nothing
    else
        return deg(b.monomials, args...)
    end
end

__combine(::Function, ::Nothing) = nothing
__combine(::Function, a) = a
__combine(f::Function, ::Nothing, args...) = __combine(f, args...)
function __combine(f::Function, a, args...)
    d = __combine(f, args...)
    if isnothing(d)
        return a
    else
        f(a, d)
    end
end
function _combine(f, args...)
    d = __combine(f, args...)
    if isnothing(d)
        error("Cannot compute `$(f)degree` as all bases are empty")
    end
    return d
end

function _combine_deg(combine, deg, d::AnyDependence, args...)
    return _combine(
        combine,
        _deg(deg, d.independent, args...),
        _deg(deg, d.dependent, args...),
    )
end

function _combine_deg(combine, deg, d::StaircaseDependence, args...)
    return _combine(
        combine,
        _deg(deg, d.standard, args...),
        _deg(deg, d.corners, args...),
        _deg(deg, d.dependent_border, args...),
        _deg(deg, d.independent_border, args...),
    )
end

function MP.mindegree(d::AbstractDependence, args...)
    return _combine_deg(min, MP.mindegree, d, args...)
end

function MP.maxdegree(d::AbstractDependence, args...)
    return _combine_deg(max, MP.maxdegree, d, args...)
end

function _ticks(d::AbstractDependence, v)
    return MP.mindegree(d, v):MP.maxdegree(d, v)
end

RecipesBase.@recipe function f(m::AnyDependence)
    vars = MP.variables(m.independent.monomials)
    t = _ticks.(Ref(m), vars)
    aspect_ratio --> :equal # defaults to `:auto`
    xticks --> t[1]
    yticks --> t[2]
    if length(t) >= 3
        zticks --> t[3]
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :circle
        label := "Independent"
        _split_exponents(m.independent)
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :rect
        label := "Dependent"
        _split_exponents(m.dependent)
    end
end

RecipesBase.@recipe function f(m::StaircaseDependence)
    vars = MP.variables(m.standard.monomials)
    t = _ticks.(Ref(m), vars)
    aspect_ratio --> :equal # defaults to `:auto`
    xticks --> t[1]
    yticks --> t[2]
    if length(t) >= 3
        zticks --> t[3]
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :circle
        label := "Standard"
        _split_exponents(m.standard)
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :rect
        label := "Corners"
        _split_exponents(m.corners)
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :rect
        label := "Dependent border"
        _split_exponents(m.dependent_border)
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :circle
        label := "Independent border"
        _split_exponents(m.independent_border)
    end
end
