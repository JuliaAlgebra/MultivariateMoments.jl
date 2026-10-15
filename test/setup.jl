using Test, LinearAlgebra, MultivariatePolynomials, MultivariateMoments
import MultivariatePolynomials as MP
import MultivariateBases as MB

if !@isdefined(Mod)
    import DynamicPolynomials
    const Mod = DynamicPolynomials
end
