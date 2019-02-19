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

function build_system(U::AbstractMatrix, mv::AbstractVector, ztol)
    # System is
    # y = [U 0] * y
    # where y = x[end:-1:1]
    # which is
    # y = U * β
    m = length(mv)
    r = size(U, 1)
    pivots = [findfirst(j -> U[i, j] != 0, 1:m) for i in 1:r]
    if any(pivots .== 0)
        keep = pivots .> 0
        pivots = pivots[keep]
        r = length(pivots)
        U = U[keep, :]
    end
    β = monovec(mv[m + 1 .- pivots]) # monovec makes sure it stays sorted, TypedPolynomials wouldn't let it sorted
    function equation(i)
        if iszero(r) # sum throws ArgumentError: reducing over an empty collection is not allowed, if r is zero
            z = zero(eltype(β)) * zero(eltype(U))
            s = z + z # For type stability
        else
            s = sum(j -> β[r+1-j] * U[j, i], 1:r)
        end
        s - mv[m+1-i]
    end
    system = filter(p -> maxdegree(p) > 0, map(equation, 1:length(mv)))
    # Type instability here :(
    if mindegree(mv) == maxdegree(mv) # Homogeneous
        projectivealgebraicset(system, Buchberger(ztol))
    else
        algebraicset(system, Buchberger(ztol))
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
    if VERSION >= v"0.7-"
        U = cholesky(M + dec.shift * I).U
    else
        U = chol(M + dec.shift * I)
    end
    σs = map(i -> (U[i, i])^2, 1:m)
    nM = maximum(σs)
    tol = nM * ranktol
    rm = findall(σs .<= tol)
    if isempty(rm)
        cM = ranktol
    else
        cM = maximum(σs[rm])
    end
    nM, cM, U[σs .> tol, :]
end
function lowrankchol(M::AbstractMatrix, dec::SVDChol, ranktol)
    if VERSION >= v"0.7-"
        F = svd(M)
    else
        F = svdfact(M)
    end
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
function computesupport!(μ::MomentMatrix, ranktol::Real, dec::LowRankChol=SVDChol())
    # We reverse the ordering so that the first columns corresponds to low order monomials
    # so that we have more chance that low order monomials are in β and then more chance
    # v[i] * β to be in μ.x
    M = getmat(μ)
    m = LinearAlgebra.checksquare(M)
    M = M[m:-1:1, m:-1:1]
    nM, cM, U = lowrankchol(M, dec, ranktol)
    W = Matrix(U)
    rref!(W, nM * cM / sqrt(m))
    #r, vals = solve_system(U', μ.x)
    μ.support = build_system(W, μ.x, sqrt(cM)) # TODO determine what is better between ranktol and sqrt(ranktol) here
end

"""
    extractatoms(ν::MomentMatrix, ranktol, [dec])

Return an `AtomicMeasure` with the atoms of `ν` if it is atomic or `nothing` if `ν` is not atomic.
The `ranktol` and `dec` parameters are passed as is to the [`lowrankchol`](@ref) function.
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
    A = Matrix{T}(undef, length(μ.x), r)
    for i in 1:r
        A[:, i] = dirac(μ.x, vars => centers[i]).a
    end
    weights = A \ μ.a
    isf = isfinite.(weights)
    weights = weights[isf]
    centers = centers[isf]
    if isempty(centers)
        nothing
    else
        AtomicMeasure(vars, WeightedDiracMeasure.(centers, weights))
    end
end
