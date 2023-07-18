module S

import MultivariatePolynomials as MP
import MultivariateBases as MB
import RecipesBase

struct MonomialDependence{B<:MB.AbstractPolynomialBasis}
    standard_monomials::B
    corners::B
    dependent_border::B
    independent_border::B
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

function _ticks(m::MonomialDependence, v)
    n = max(
        MP.maxdegree(m.standard_monomials.monomials, v),
        MP.maxdegree(m.corners.monomials, v),
        MP.maxdegree(m.dependent_border.monomials, v),
        MP.maxdegree(m.independent_border.monomials, v),
    )
    return 0:n
end

RecipesBase.@recipe function f(m::MonomialDependence)
    vars = MP.variables(m.standard_monomials.monomials)
    t = _ticks.(Ref(m), vars)
    aspect_ratio --> :equal # defaults to :auto
    xticks --> t[1]
    yticks --> t[2]
    if length(t) >= 3
        zticks --> t[3]
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :circle
        label := "Standard"
        _split_exponents(m.standard_monomials)
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
        label := "Dependent"
        _split_exponents(m.dependent_border)
    end
    RecipesBase.@series begin
        seriestype --> :scatter
        markershape --> :circle
        label := "Independent"
        _split_exponents(m.independent_border)
    end
end

end
