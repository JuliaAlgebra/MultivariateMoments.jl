# The content of this file is heavily inspired from MultivariateSeries.jl

struct ZeroDimensionalVariety{T} <: AbstractAlgebraicSet
    elements::Vector{Vector{T}}
end
SemialgebraicSets.is_zero_dimensional(::ZeroDimensionalVariety) = true
Base.length(v::ZeroDimensionalVariety) = length(v.elements)
Base.iterate(v::ZeroDimensionalVariety, args...) = iterate(v.elements, args...)

# Decomposition of the pencil of matrices
function decompose(
    H::Vector{Matrix{T}},
    λ::Vector,
    rank_check::RankCheck,
    multiplication_matrices_solver::SemialgebraicSets.AbstractMultiplicationMatricesSolver,
) where {T}
    H0 = sum(H[i] * λ[i] for i in eachindex(λ))

    # H0 = U * Diag(σ) * V'
    U, σ, V = svd(H0)
    r = rank_from_singular_values(σ, rank_check)

    Σi = Diagonal(inv.(σ[1:r]))

    M = Matrix{T}[Σi * (U[:, 1:r]') * H[i] * (V[:, 1:r]) for i in 2:length(H)]
    mult = SemialgebraicSets.MultiplicationMatrices(M)

    if r > 1
        return SemialgebraicSets.solve(mult, multiplication_matrices_solver)
    else
        return [[M[i][end, end] for i in eachindex(M)]]
    end
end

struct FlatExtension{
    MMS<:SemialgebraicSets.AbstractMultiplicationMatricesSolver,
}
    multiplication_matrices_solver::MMS
end
function FlatExtension()
    return FlatExtension(ReorderedSchurMultiplicationMatricesSolver{Float64}())
end

function hankel(μ::Measure{T}, rows, cols) where {T}
    return T[
        moment_value(μ, rows[i] * cols[j]) for i in eachindex(rows),
        j in eachindex(cols)
    ]
end

function compute_support!(
    ν::MomentMatrix{T},
    rank_check::RankCheck,
    solver::FlatExtension,
) where {T}
    μ = measure(ν)
    d = MP.maxdegree(μ.x)
    v = MP.variables(μ)
    d0 = div(d - 1, 2)
    d1 = d - 1 - d0
    B0 = MP.monomials(v, 0:d0)
    B1 = MP.monomials(v, 0:d1)

    H = Matrix{T}[hankel(μ, B0, B1)]
    for x in v
        push!(H, hankel(μ, B0, [b * x for b in B1]))
    end
    λ = [one(T)]
    ν.support = ZeroDimensionalVariety(
        decompose(H, λ, rank_check, solver.multiplication_matrices_solver),
    )
    return
end
