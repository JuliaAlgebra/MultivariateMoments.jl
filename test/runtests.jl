using MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
using MultivariateMoments
using SemialgebraicSets

using Test
using LinearAlgebra

include("symmatrix.jl")
include("hermitian_matrix.jl")

# Taken from JuMP/test/solvers.jl
function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

if try_import(:DynamicPolynomials)
    Mod = DynamicPolynomials
    include("commutativetests.jl")
    #include("noncommutativetests.jl")
end

if try_import(:TypedPolynomials)
    Mod = TypedPolynomials
    include("commutativetests.jl")
end
