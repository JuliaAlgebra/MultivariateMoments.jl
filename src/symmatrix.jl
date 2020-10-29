export SymMatrix

"""
    struct SymMatrix{T} <: AbstractMatrix{T}
        Q::Vector{T}
        n::Int
    end

Symmetric ``n \\times n`` matrix storing the vectorized upper triangular
part of the matrix in the `Q` vector (this corresponds to the vectorized
form of `MathOptInterface.PositiveSemidefiniteConeTriangle`).
It implement the `AbstractMatrix` interface except for `setindex!` as it might
break its symmetry. The [`symmetric_setindex!`](@ref) function should be used
instead.
"""
struct SymMatrix{T} <: AbstractMatrix{T}
    Q::Vector{T}
    n::Int
end

_undef_sym(T::Type, n) = Vector{T}(undef, trimap(n, n))
function SymMatrix{T}(::UndefInitializer, dims::Dims{2}) where T
    dims[1] != dims[2] && error("Expected same dimension for `SymMatrix`, got `$(dims)`.")
    n = dims[1]
    return SymMatrix(_undef_sym(T, n), n)
end
Base.similar(Q::SymMatrix, T::Type, dims::Dims{2}) = similar(SymMatrix{T}, dims)
Base.copy(Q::SymMatrix) = SymMatrix(copy(Q.Q), Q.n)
Base.map(f::Function, Q::SymMatrix) = SymMatrix(map(f, Q.Q), Q.n)

# i <= j
#function trimap(i, j, n) # It was used when Q was the lower triangular part
#    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
#end

# i <= j
trimap(i, j) = div(j * (j - 1), 2) + i

function trimat(::Type{T}, f, n, σ) where {T}
    Q = _undef_sym(T, n)
    for j in 1:n
        for i in 1:j
            Q[trimap(i, j)] = f(σ[i], σ[j])
        end
    end
    return SymMatrix(Q, n)
end

Base.size(Q::SymMatrix) = (Q.n, Q.n)

"""
    symmetric_setindex!(Q::SymMatrix, value, i::Integer, j::Integer)

Set `Q[i, j]` and `Q[j, i]` to the value `value`.
"""
function symmetric_setindex!(Q::SymMatrix, value, i::Integer, j::Integer)
    Q.Q[trimap(min(i, j), max(i, j))] = value
end

function Base.getindex(Q::SymMatrix, i::Integer, j::Integer)
    return Q.Q[trimap(min(i, j), max(i, j))]
end
Base.getindex(Q::SymMatrix, I::Tuple) = Q[I...]
Base.getindex(Q::SymMatrix, I::CartesianIndex) = Q[I.I]
