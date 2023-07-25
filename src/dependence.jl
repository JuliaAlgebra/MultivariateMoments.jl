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
    if !isempty(d.independent.monomials)
        println(io, "with independent basis:")
        println(io, d.independent)
    end
    if !isempty(d.dependent.monomials)
        println(io, "and dependent basis:")
        println(io, d.dependent)
    end
    return
end

"""
    struct StaircaseDependence{B<:MB.AbstractPolynomialBasis} <: AbstractDependence
        trivial_standard::B
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
    trivial_standard::B
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
    if !isempty(d.trivial_standard.monomials)
        println(io, "with trivial standard basis:")
        println(io, d.trivial_standard)
    end
    if !isempty(d.standard.monomials)
        println(io, "with standard basis:")
        println(io, d.standard)
    end
    if !isempty(d.corners.monomials)
        println(io, "and corners basis:")
        println(io, d.corners)
    end
    if !isempty(d.dependent_border.monomials)
        println(io, "and dependent border basis:")
        println(io, d.dependent_border)
    end
    if !isempty(d.independent_border.monomials)
        println(io, "and independent border basis:")
        print(io, d.independent_border)
    end
    return
end

function Base.convert(::Type{AnyDependence}, d::StaircaseDependence)
    dependent = MB.merge_bases(d.corners, d.dependent_border)[1]
    return AnyDependence(d.standard, dependent)
end

function AnyDependence(
    is_dependent::Function,
    basis::MB.MonomialBasis{M},
) where {M}
    independent = M[]
    dependent = M[]
    for i in eachindex(basis.monomials)
        mono = basis.monomials[i]
        if is_dependent(i)
            push!(dependent, mono)
        else
            push!(independent, mono)
        end
    end
    return AnyDependence(
        MB.MonomialBasis(independent),
        MB.MonomialBasis(dependent),
    )
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
    trivial_standard = M[]
    standard = M[]
    corners = M[]
    vars = MP.variables(basis)
    for mono in MP.monomials(vars, 0:MP.maxdegree(basis.monomials))
        # This sieve of [LLR08, Algorithm 1] is a performance improvement but not only.
        # It also ensures that the standard monomials have the "staircase structure".
        if !any(Base.Fix2(MP.divides, mono), corners)
            i = _index(basis, mono)
            if isnothing(i)
                push!(trivial_standard, mono)
            elseif !is_dependent(i)
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
                i = _index(basis, border)
                if isnothing(i) || !is_dependent(i)
                    push!(independent_border, border)
                else
                    push!(dependent_border, border)
                end
            end
        end
    end
    return StaircaseDependence(
        MB.MonomialBasis(trivial_standard),
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

__reduce(::Function, ::Nothing) = nothing
__reduce(::Function, a) = a
__reduce(f::Function, ::Nothing, args...) = __reduce(f, args...)
function __reduce(f::Function, a, args...)
    d = __reduce(f, args...)
    if isnothing(d)
        return a
    else
        f(a, d)
    end
end
function _reduce(f, args...)
    d = __reduce(f, args...)
    if isnothing(d)
        error("Cannot compute `$(f)degree` as all bases are empty")
    end
    return d
end

function _map_reduce(combine, deg, d::AnyDependence, args...)
    return _reduce(
        combine,
        _deg(deg, d.independent, args...),
        _deg(deg, d.dependent, args...),
    )
end

function _map_reduce(combine, deg, d::StaircaseDependence, args...)
    return _reduce(
        combine,
        _deg(deg, d.trivial_standard, args...),
        _deg(deg, d.standard, args...),
        _deg(deg, d.corners, args...),
        _deg(deg, d.dependent_border, args...),
        _deg(deg, d.independent_border, args...),
    )
end

function _reduce_variables(v, w)
    return MP.variables(prod(v) * prod(w))
end

function MP.variables(d::AbstractDependence)
    return _map_reduce(_reduce_variables, MP.variables, d)
end

function MP.mindegree(d::AbstractDependence, args...)
    return _map_reduce(min, MP.mindegree, d, args...)
end

function MP.maxdegree(d::AbstractDependence, args...)
    return _map_reduce(max, MP.maxdegree, d, args...)
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
    if !isempty(m.independent.monomials)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> :circle
            label := "Independent"
            _split_exponents(m.independent)
        end
    end
    if !isempty(m.dependent.monomials)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> :rect
            label := "Dependent"
            _split_exponents(m.dependent)
        end
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
    if !isempty(m.trivial_standard.monomials)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> :circle
            label := "Trivial standard"
            _split_exponents(m.trivial_standard)
        end
    end
    if !isempty(m.standard.monomials)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> :circle
            label := "Standard"
            _split_exponents(m.standard)
        end
    end
    if !isempty(m.corners.monomials)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> :rect
            label := "Corners"
            _split_exponents(m.corners)
        end
    end
    if !isempty(m.dependent_border.monomials)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> :rect
            label := "Dependent border"
            _split_exponents(m.dependent_border)
        end
    end
    if !isempty(m.independent_border.monomials)
        RecipesBase.@series begin
            seriestype --> :scatter
            markershape --> :circle
            label := "Independent border"
            _split_exponents(m.independent_border)
        end
    end
end
