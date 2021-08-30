export extractatoms
export LowRankChol, ShiftChol, SVDChol

using RowEchelon
using SemialgebraicSets

# Solve the system
# y = [U 0] * y
# where y = x[end:-1:1]
# which is
# y = U * β
# The code for solving the system with reordered schur has been moved to SemialgebraicSets.jl
#function solve_system(U, x)
#    m, r = size(U)
#    @assert m == length(x)
#    n = nvariables(x)
#    v = variables(x)
#    pivots = [findfirst(j -> U[j, i] != 0, 1:m) for i in 1:r]
#    if any(pivots .== 0)
#        keep = pivots .> 0
#        pivots = pivots[keep]
#        r = length(pivots)
#        U = U[:, keep]
#    end
#    β = x[m+1-pivots]
#    function multisearch_check(y)
#        idxs = multisearch(x, y)
#        if any(idxs .== 0)
#            error("Missing monomials $(y[idxs .== 0]) in $(x)")
#        end
#        idxs
#    end
#    Ns = [U[m+1-reverse(multisearch_check(v[i] * β)), :] for i in 1:n]
#    λ = rand(n)
#    λ /= sum(λ)
#    N = sum(λ .* Ns)
#    Z = schurfact(N)[:Z]
#    vals = [Vector{Float64}(n) for j in 1:r]
#    for j in 1:r
#        qj = Z[:, j]
#        for i in 1:n
#            vals[j][i] = dot(qj, Ns[i] * qj)
#        end
#    end
#    r, vals
#end

function build_system(U::AbstractMatrix, basis::MB.MonomialBasis, ztol, args...)
    # System is
    # y = [U 0] * y
    # where y = x[end:-1:1]
    # which is
    # y = U * β
    m = length(basis)
    r = size(U, 1)
    pivots = [findfirst(j -> U[i, j] != 0, 1:m) for i in 1:r]
    if any(pivots .== 0)
        keep = pivots .> 0
        pivots = pivots[keep]
        r = length(pivots)
        U = U[keep, :]
    end
    monos = basis.elements
    β = monovec(monos[m + 1 .- pivots]) # monovec makes sure it stays sorted, TypedPolynomials wouldn't let it sorted
    function equation(i)
        if iszero(r) # sum throws ArgumentError: reducing over an empty collection is not allowed, if r is zero
            z = zero(eltype(β)) * zero(eltype(U))
            s = z + z # For type stability
        else
            s = sum(j -> β[r+1-j] * U[j, i], 1:r)
        end
        s - monos[m+1-i]
    end
    system = filter(p -> maxdegree(p) > 0, map(equation, 1:length(monos)))
    # Type instability here :(
    if mindegree(monos) == maxdegree(monos) # Homogeneous
        projectivealgebraicset(system, Buchberger(ztol), args...)
    else
        algebraicset(system, Buchberger(ztol), args...)
    end
end

"""
    LowRankChol

Method for computing a ``r \\times n`` matrix `U` of a ``n \\times n`` rank ``r`` psd matrix `Q` such that `Q = U'U`.
"""
abstract type LowRankChol end

"""
    ShiftChol <: LowRankChol

Shift the matrix by `shift` times the identity matrix before cholesky.
"""
struct ShiftChol{T} <: LowRankChol
    shift::T
end

"""
    SVDChol <: LowRankChol

Use SVD decomposition.
"""
struct SVDChol <: LowRankChol end

"""
    MultivariateMoments.lowrankchol(Q::AbstractMatrix, dec::LowRankChol, ranktol)

Returns a ``r \\times n`` matrix ``U`` of a ``n \\times n`` rank ``r`` positive semidefinite matrix ``Q`` such that ``Q = U^\\top U``.
The rank of ``Q`` is the number of singular values larger than `ranktol```{} \\cdot \\sigma_1`` where ``\\sigma_1`` is the largest singular value.
"""
function lowrankchol end

function lowrankchol(M::AbstractMatrix, dec::ShiftChol, ranktol)
    m = LinearAlgebra.checksquare(M)
    U = cholesky(M + dec.shift * I).U
    σs = map(i -> (U[i, i])^2, 1:m)
    nM = maximum(σs)
    tol = nM * ranktol
    rm = findall(σs .<= tol)
    if isempty(rm)
        cM = ranktol
    else
        cM = maximum(σs[rm]) / nM
    end
    nM, cM, U[σs .> tol, :]
end
function lowrankchol(M::AbstractMatrix, dec::SVDChol, ranktol)
    F = svd(M)
    nM = F.S[1] # norm of M
    tol = nM * ranktol
    r = something(findfirst(σ2 -> σ2 <= tol, F.S), 0)
    if iszero(r)
        cM = ranktol
        r = length(F.S)
    else
        cM = F.S[r] / nM
        r -= 1
    end
    nM, cM, (F.U[:, 1:r] * Diagonal(sqrt.(F.S[1:r])))'
end

"""
    MultivariateMoments.computesupport!(ν::MomentMatrix, ranktol, [dec])

Computes the `support` field of `ν`.
The `ranktol` and `dec` parameters are passed as is to the [`lowrankchol`](@ref) function.
"""
function computesupport!(μ::MomentMatrix, ranktol::Real, dec::LowRankChol, args...)
    # We reverse the ordering so that the first columns corresponds to low order monomials
    # so that we have more chance that low order monomials are in β and then more chance
    # v[i] * β to be in μ.x
    M = getmat(μ)
    m = LinearAlgebra.checksquare(M)
    M = M[m:-1:1, m:-1:1]
    nM, cM, U = lowrankchol(M, dec, ranktol)
    W = Matrix(U)
    # If M is multiplied by λ, W is multiplied by √λ
    # so we take √||M|| = √nM
    rref!(W, √(nM) * cM / sqrt(m))
    #r, vals = solve_system(U', μ.x)
    μ.support = build_system(W, μ.basis, √cM, args...) # TODO determine what is better between ranktol and sqrt(ranktol) here
end

function computesupport!(μ::MomentMatrix, ranktol::Real, args...)
    return computesupport!(μ::MomentMatrix, ranktol::Real, SVDChol(), args...)
end

"""
    extractatoms(ν::MomentMatrix, ranktol, [dec::LowRankChol], [solver::SemialgebraicSets.AbstractAlgebraicSolver])

Return an `AtomicMeasure` with the atoms of `ν` if it is atomic or `nothing` if
`ν` is not atomic. The `ranktol` and `dec` parameters are passed as is to the
[`lowrankchol`](@ref) function. By default, `dec` is an instance of
[`SVDChol`](@ref). The extraction relies on the solution of a system of
algebraic equations. using `solver`. For instance, given a
[`MomentMatrix`](@ref), `μ`, the following extract atoms using a `ranktol` of
`1e-4` for the low-rank decomposition and homotopy continuation to solve the
obtained system of algebraic equations.
```julia
using HomotopyContinuation
solver = SemialgebraicSetsHCSolver(; compile = false)
atoms = extractatoms(ν, 1e-4, solver)
```
If no solver is given, the default solver of SemialgebraicSets is used which
currently computes the Gröbner basis, then the multiplication matrices and
then the Schur decomposition of a random combination of these matrices.
For floating point arithmetics, homotopy continuation is recommended as it is
more numerically stable than Gröbner basis computation.
"""
function extractatoms(ν::MomentMatrix{T}, ranktol, args...) where T
    computesupport!(ν, ranktol, args...)
    supp = ν.support
    if !iszerodimensional(supp)
        return nothing
    end
    centers = collect(supp)
    r = length(centers)
    # Determine weights
    μ = measure(ν)
    vars = variables(μ)
    A = Matrix{T}(undef, length(μ.basis), r)
    for i in 1:r
        A[:, i] = dirac(μ.basis.elements, vars => centers[i]).values
    end
    weights = A \ μ.values
    isf = isfinite.(weights)
    weights = weights[isf]
    centers = centers[isf]
    if isempty(centers)
        nothing
    else
        AtomicMeasure(vars, WeightedDiracMeasure.(centers, weights))
    end
end
