export SymMatrix, MatMeasure, getmat, matmeasure, AtomicMeasure, extractatoms
using RowEchelon
using SemialgebraicSets

struct SymMatrix{T} <: AbstractMatrix{T}
    Q::Vector{T}
    n
end

# i < j
function trimap(i, j, n)
    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
end

function trimat(::Type{T}, f, n, σ) where {T}
    Q = Vector{T}(trimap(n, n, n))
    for i in 1:n
        for j in i:n
            Q[trimap(i, j, n)] = f(σ[i], σ[j])
        end
    end
    SymMatrix{T}(Q, n)
end

Base.size(Q::SymMatrix) = (Q.n, Q.n)

function Base.getindex(Q::SymMatrix, i, j)
    Q.Q[trimap(min(i, j), max(i, j), Q.n)]
end
function Base.getindex(Q::SymMatrix, k)
    i, j = divrem(k-1, Q.n)
    Q[i+1, j+1]
end
Base.getindex(Q::SymMatrix, I::Tuple) = Q[I...]
Base.getindex(Q::SymMatrix, I::CartesianIndex) = Q[I.I]

mutable struct MatMeasure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractMeasureLike{T}
    Q::SymMatrix{T}
    x::MVT
    support::Nullable{AlgebraicSet}
end
MatMeasure{T, MT, MVT}(Q::SymMatrix{T}, x::MVT) where {T, MT, MVT} = MatMeasure{T, MT, MVT}(Q, x, nothing)

MP.variables(μ::MatMeasure) = variables(μ.x)
MP.nvariables(μ::MatMeasure) = nvariables(μ.x)

function MatMeasure{T}(f::Function, x::AbstractVector{MT}, σ) where {T, MT<:AbstractMonomial}
    MatMeasure{T, MT, monovectype(x)}(trimat(T, f, length(x), σ), x)
end
MatMeasure{T}(f::Function, x::AbstractVector, σ) where T = MatMeasure{T}(f, monomial.(x), σ)
function MatMeasure{T}(f::Function, x::AbstractVector) where T
    σ, X = sortmonovec(x)
    MatMeasure{T}(f, X, σ)
end
MatMeasure(f::Function, x) = MatMeasure{Base.promote_op(f, Int, Int)}(f, x)

function matmeasure(μ::Measure{T}, X) where T
    function getmom(i, j)
        x = X[i] * X[j]
        for m in moments(μ)
            if monomial(m) == x
                return value(m)
            end
        end
        throw(ArgumentError("μ does not have the moment $(x)"))
    end
    MatMeasure{T}(getmom, X)
end

function MatMeasure(Q::AbstractMatrix{T}, x, σ) where T
    MatMeasure{T}((i,j) -> Q[σ[i], σ[j]], x)
end
function MatMeasure(Q::AbstractMatrix, x)
    σ, X = sortmonovec(x)
    MatMeasure(Q, X, σ)
end

function getmat{C, T}(μ::MatMeasure{C, T})
    Matrix(μ.Q)
end

type AtomicMeasure{T, V}
    v::V # Vector/Tuple of variables
    λ::Vector{T} # The measure is sum λ_i * δ_{support_i}
    support::Vector{Vector{T}} # Elements of the finite support
end
function AtomicMeasure{V, S, T}(vars::V, λ::Vector{S}, support::Vector{Vector{T}})
    AtomicMeasure{promote_type(S, T), V}(vars, λ, support)
end

function Measure(μ::AtomicMeasure{T}, x::AbstractVector{TT}) where {T, TT}
    Measure{T, monomialtype(TT), monovectype(x)}(μ, x)
end
function Measure{T, MT, MVT}(μ::AtomicMeasure{T}, x) where {T, MT, MVT}
    X = monovec(x)
    sum(μ.λ[i] * dirac(X, μ.v=>μ.support[i]) for i in 1:length(μ.λ))
end

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

function build_system(U::AbstractMatrix, mv::AbstractVector)
    # System is
    # y = [U 0] * y
    # where y = x[end:-1:1]
    # which is
    # y = U * β
    m = length(mv)
    equation(i) = sum(j -> mv[m+1-j] * U[j, i], 1:size(U, 1)) - mv[m+1-i]
    system = filter(p -> maxdegree(p) > 0, map(equation, 1:length(mv)))
    algebraicset(system)
end

function computesupport!(μ::MatMeasure, tol::Real, shift::Real, ɛ::Real=-1)
    # We reverse the ordering so that the first columns corresponds to low order monomials
    # so that we have more chance that low order monomials are in β and then more chance
    # v[i] * β to be in μ.x
    M = getmat(μ)[end:-1:1, end:-1:1]
    m = size(M, 1)
    #r = rank(M, tol)
    U = chol(M + shift * eye(m))
    V = U[map(i -> abs(U[i, i]) > tol, 1:m), :]
    r = size(V, 1)
#   F = svdfact(M)
#   S = F.S
#   r = sum(F.S .> tol)
#   V = F.U[:, 1:r] .* repmat(sqrt.(S[1:r])', size(F.U, 1), 1)
    rref!(V, ɛ == -1 ? sqrt(eps(norm(V, Inf))) : ɛ)
    #r, vals = solve_system(V', μ.x)
    μ.support = build_system(V, μ.x)
end

function extractatoms(μ::MatMeasure, tol::Real, shift::Real, ɛ::Real=-1)
    M = getmat(μ)[end:-1:1, end:-1:1]
    computesupport!(μ, tol, shift, ɛ)
    supp = get(μ.support)
    if !iszerodimensional(supp)
        error("Cannot extract atoms of Measure with non zero-dimensional support")
    end
    vals = collect(supp)
    r = length(vals)
    # Determine weights
    Ms = similar(M, r, r)
    v = variables(μ)
    for i in 1:r
        vi = dirac(μ.x, v => vals[i])
        Ms[:, i] = vi.a[end:-1:end-r+1] * vi.a[end]
    end
    λ = Ms \ M[1:r, 1]
    AtomicMeasure(v, λ, vals)
end

function permcomp(f, m)
    picked = IntSet()
    for i in 1:m
        k = 0
        for j in 1:m
            if !(j in picked) && f(i, j)
                k = j
                break
            end
        end
        if k == 0
            return false
        end
        push!(picked, k)
    end
    true
end
function Base.isapprox(μ::AtomicMeasure, ν::AtomicMeasure; kws...)
    m = length(μ.λ)
    if length(ν.λ) != m
        false
    else
        permcomp((i, j) -> isapprox(μ.λ[i], ν.λ[j]; kws...) && isapprox(μ.support[i], ν.support[j]; kws...), m)
    end
end
