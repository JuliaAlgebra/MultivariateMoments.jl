export FlatExtension, IterativeDiagonalization

struct ZeroDimensionalVariety{T} <: AbstractAlgebraicSet
    elements::Vector{Vector{T}}
end
SemialgebraicSets.iszerodimensional(::ZeroDimensionalVariety) = true
Base.length(v::ZeroDimensionalVariety) = length(v.elements)
Base.iterate(v::ZeroDimensionalVariety, args...) = iterate(v.elements, args...)

# norm of off diagonal terms of a square matrix
function norm_off(M)
    if size(M[1], 1) > 1
        return sqrt(sum(abs2(M[i,j]) + abs2(M[j,i]) for i in 1:size(M,1) for j in i+1:size(M,1)))
    else
        return 0.0
    end
end

function diagonalization_iter(D)
    n = size(D[1],1)
    s = length(D)
    
    X = fill(zero(D[1][1,1]),n,n)
    Y = fill(zero(D[1][1,1]),n,n)

    A = fill(zero(D[1][1,1]),s,2)
    b = fill(zero(D[1][1,1]),s)
    for i in 1:n
        for j in 1:n
            if i != j
                for k in 1:s
                    A[k,1] = D[k][i,i]
                    A[k,2] = D[k][j,j]
                    b[k]   = -D[k][i,j]
                end
                v = A\b
                X[i,j] =  v[1]
                Y[i,j] =  v[2]
            end
        end
    end
    for i in 1:n
        X[i,i]=1
        Y[i,i]=1
    end
    return X, Y
end

struct IterativeDiagonalization{T} <: SemialgebraicSets.AbstractMultiplicationMatricesSolver
    max_iter::Int
    ε::T
    tol::T
end
# These were the values in MultivariateSeries/diagonalization.jl
IterativeDiagonalization() = IterativeDiagonalization(10, 1e-3, 5e-2)

function SemialgebraicSets.solvemultiplicationmatrices(M::Vector{Matrix{C}}, solver::IterativeDiagonalization) where C
    n  = length(M)
    r  = size(M[1], 1)

    M1 = sum(M[i]*randn(Float64) for i in 1:n)
    E  = eigvecs(M1)

    F  = inv(E)
    
    D  = vcat([Matrix{C}(I,r,r)], [F*M[i]*E for i in 1:length(M)])
    err = sum(norm_off.(D))
    Δ = sum(norm.(D))

    nit = 0

    if err / Δ > solver.tol
        Δ = err
        while nit < solver.max_iter && Δ > solver.ε
            err0 = err
            X, Y = diagonalization_iter(D)
            D = [Y*D[i]*X for i in 1:length(D)]
            E = E * X
            F = Y * F
            nit += 1
            err = sum(norm_off.(D))
            Δ = err0 - err
        end
    end
    
    return [[D[j+1][i,i] / D[1][i,i] for j in 1:n] for i in 1:r]
end

# Inspired from MultivariateSeries.jl
# Decomposition of the pencil of matrices
function decompose(
    H::Vector{Matrix{T}},
    λ::Vector,
    rank_check::RankCheck,
    multiplication_matrices_solver::SemialgebraicSets.AbstractMultiplicationMatricesSolver,
) where T
    H0 = sum(H[i] * λ[i] for i in eachindex(λ))

    # H0 = U * Diag(σ) * V'
    U, σ, V = svd(H0)
    r = rank_from_singular_values(σ, rank_check)

    Σi = Diagonal(inv.(σ[1:r]))

    M = Matrix{T}[Σi * (U[:,1:r]') * H[i] * (V[:,1:r]) for i in 2:length(H)]

    if r > 1
        return SemialgebraicSets.solvemultiplicationmatrices(M, multiplication_matrices_solver)
    else
        return [[M[i][end, end] for i in eachindex(M)]]
    end
end

struct FlatExtension{MMS<:SemialgebraicSets.AbstractMultiplicationMatricesSolver}
    multiplication_matrices_solver::MMS
end
FlatExtension() = FlatExtension(ReorderedSchurMultiplicationMatricesSolver{Float64}())

function hankel(μ::Measure{T}, rows, cols) where T
    return T[moment_value(μ, rows[i] * cols[j]) for i in eachindex(rows), j in eachindex(cols)]
end

function computesupport!(ν::MomentMatrix{T}, rank_check::RankCheck, solver::FlatExtension) where T
    μ = measure(ν)
    d = maxdegree(μ.x)
    v = variables(μ)
    d0 = div(d - 1, 2)
    d1 = d - 1 - d0
    B0 = monomials(v, 0:d0)
    B1 = monomials(v, 0:d1)

    H = Matrix{T}[hankel(μ, B0, B1)]
    for x in v
        push!(H, hankel(μ, B0, [b*x for b in B1]))
    end
    λ = [one(T)]
    ν.support = ZeroDimensionalVariety(decompose(H, λ, rank_check, solver.multiplication_matrices_solver))
    return
end
