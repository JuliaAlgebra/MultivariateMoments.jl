export atomic_measure
export MomentMatrixWeightSolver, MomentVectorWeightSolver

using RowEchelon
using SemialgebraicSets

"""
    struct MacaulayNullspace{T,MT<:AbstractMatrix{T},BT}
        matrix::MT
        basis::BT
    end

Right nullspace with rows indexed by `basis` of a Macaulay matrix with columns
indexed by `basis`.
The value of `matrix[i, j]` should be interpreted as the value of the `i`th
element of `basis` for the `j`th generator of the nullspace.
"""
struct MacaulayNullspace{T,MT<:AbstractMatrix{T},BT}
    matrix::MT
    basis::BT
    accuracy::T
end

"""
    MultivariateMoments.compute_support!(ν::MomentMatrix, rank_check, [dec])

Computes the `support` field of `ν`.
The `rank_check` and `dec` parameters are passed as is to the [`low_rank_cholesky`](@ref) function.
"""
function compute_support!(μ::MomentMatrix, rank_check::RankCheck, dec::LowRankLDLTAlgorithm, nullspace_solver=Echelon(), args...)
    # Ideally, we should determine the pivots with a greedy sieve algorithm [LLR08, Algorithm 1]
    # so that we have more chance that low order monomials are in β and then more chance
    # so that the pivots form an order ideal. We just use `rref` which does not implement the sieve
    # v[i] * β to be in μ.x
    # but maybe it's sufficient due to the shift structure of the matrix.
    #
    # [LLR08] Lasserre, Jean Bernard and Monique Laurent, and Philipp Rostalski.
    # "Semidefinite characterization and computation of zero-dimensional real radical ideals."
    # Foundations of Computational Mathematics 8 (2008): 607-647.
    M = value_matrix(μ)
    chol = low_rank_cholesky(M, dec, rank_check)
    @assert size(chol.L, 1) == LinearAlgebra.checksquare(M)
    μ.support = solve(MacaulayNullspace(chol.L, μ.basis, recommended_rtol(chol)), nullspace_solver, args...)
end

function compute_support!(μ::MomentMatrix, rank_check::RankCheck, args...)
    return compute_support!(μ::MomentMatrix, rank_check::RankCheck, SVDCholesky(), args...)
end

# Determines weight

"""
    struct MomentMatrixWeightSolver
        rtol::T
        atol::T
    end

Given a moment matrix `ν` and the atom centers,
determine the weights by solving a linear system over all the moments
of the moment matrix, keeping duplicates (e.g., entries corresponding to the same monomial).

If the moment values corresponding to the same monomials are known to be equal
prefer [`MomentVectorWeightSolver`](@ref) instead.
"""
struct MomentMatrixWeightSolver
end

function solve_weight(ν::MomentMatrix{T}, centers, ::MomentMatrixWeightSolver) where {T}
    vars = MP.variables(ν)
    A = Matrix{T}(undef, length(ν.Q.Q), length(centers))
    vbasis = vectorized_basis(ν)
    for i in eachindex(centers)
        η = dirac(vbasis.monomials, vars => centers[i])
        A[:, i] = moment_matrix(η, ν.basis.monomials).Q.Q
    end
    return A \ ν.Q.Q
end

"""
    struct MomentVectorWeightSolver{T}
        rtol::T
        atol::T
    end

Given a moment matrix `ν` and the atom centers, first convert the moment matrix
to a vector of moments, using [`measure(ν; rtol=rtol, atol=atol)`](@ref measure)
and then determine the weights by solving a linear system over the monomials obtained.

If the moment values corresponding to the same monomials can have small differences,
[`measure`](@ref) can throw an error if `rtol` and `atol` are not small enough.
Alternatively to tuning these tolerances [`MomentVectorWeightSolver`](@ref) can be used instead.
"""
struct MomentVectorWeightSolver{T}
    rtol::T
    atol::T
end
function MomentVectorWeightSolver{T}(; rtol=Base.rtoldefault(T), atol=zero(T)) where {T}
    return MomentVectorWeightSolver{T}(rtol, atol)
end
function MomentVectorWeightSolver(; rtol=nothing, atol=nothing)
    if rtol === nothing && atol === nothing
        return MomentVectorWeightSolver{Float64}()
    elseif rtol !== nothing
        if atol === nothing
            return MomentVectorWeightSolver{typeof(rtol)}(; rtol=rtol)
        else
            return MomentVectorWeightSolver{typeof(rtol)}(; rtol=rtol, atol=atol)
        end
    else
        return MomentVectorWeightSolver{typeof(atol)}(; atol=atol)
    end
end

function solve_weight(ν::MomentMatrix{T}, centers, solver::MomentVectorWeightSolver) where {T}
    μ = measure(ν; rtol=solver.rtol, atol=solver.atol)
    vars = MP.variables(μ)
    A = Matrix{T}(undef, length(μ.x), length(centers))
    for i in eachindex(centers)
        A[:, i] = dirac(μ.x, vars => centers[i]).a
    end
    return A \ μ.a
end

"""
    atomic_measure(ν::MomentMatrix, rank_check::RankCheck, [dec::LowRankLDLTAlgorithm], [solver::SemialgebraicSets.AbstractAlgebraicSolver])

Return an `AtomicMeasure` with the atoms of `ν` if it is atomic or `nothing` if
`ν` is not atomic. The `rank_check` and `dec` parameters are passed as is to the
[`low_rank_cholesky`](@ref) function. By default, `dec` is an instance of
[`SVDCholesky`](@ref). The extraction relies on the solution of a system of
algebraic equations. using `solver`. For instance, given a
[`MomentMatrix`](@ref), `μ`, the following extract atoms using a `rank_check` of
`1e-4` for the low-rank decomposition and homotopy continuation to solve the
obtained system of algebraic equations.
```julia
using HomotopyContinuation
solver = SemialgebraicSetsHCSolver(; compile = false)
atoms = atomic_measure(ν, 1e-4, solver)
```
If no solver is given, the default solver of SemialgebraicSets is used which
currently computes the Gröbner basis, then the multiplication matrices and
then the Schur decomposition of a random combination of these matrices.
For floating point arithmetics, homotopy continuation is recommended as it is
more numerically stable than Gröbner basis computation.
"""
function atomic_measure(ν::MomentMatrix, rank_check::RankCheck, args...; weight_solver = MomentMatrixWeightSolver())
    compute_support!(ν, rank_check, args...)
    supp = ν.support
    if isnothing(supp)
        return
    end
    if !is_zero_dimensional(supp)
        return
    end
    centers = collect(supp)
    weights = solve_weight(ν, centers, weight_solver)
    isf = isfinite.(weights)
    weights = weights[isf]
    centers = centers[isf]
    if isempty(centers)
        nothing
    else
        AtomicMeasure(MP.variables(ν), WeightedDiracMeasure.(centers, weights))
    end
end

function atomic_measure(ν::MomentMatrix, ranktol, args...; kws...)
    return atomic_measure(ν, LeadingRelativeRankTol(ranktol), args...; kws...)
end
