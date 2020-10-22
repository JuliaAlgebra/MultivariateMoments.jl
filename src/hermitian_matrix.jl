export VectorizedHermitianMatrix

"""
    struct VectorizedHermitianMatrix{T} <: AbstractMatrix{T}
        Q::Vector{T}
        n::Int
    end

Hermitian ``n \\times n`` matrix storing the vectorized upper triangular
real part of the matrix followed by the vectorized upper triangular
imaginary part in the `Q` vector (this corresponds to the vectorized
form of `ComplexOptInterface.HermitianPositiveSemidefiniteConeTriangle`).
It implement the `AbstractMatrix` interface except for `setindex!` as it might
break its symmetry. The [`symmetric_setindex!`](@ref) function should be used
instead.
"""
struct VectorizedHermitianMatrix{T} <: AbstractMatrix{T}
    Q::Vector{T}
    n::Int
end

Base.copy(Q::VectorizedHermitianMatrix) = VectorizedHermitianMatrix(copy(Q.Q), Q.n)
function Base.map(f::Function, Q::VectorizedHermitianMatrix)
    if Q.n <= 1
        return VectorizedHermitianMatrix(map(f, Q.Q), Q.n)
    else
        el = f(Q[1, 2])
        ret = VectorizedHermitianMatrix(similar(Q.Q, typeof(real(el))), Q.n)
        for j in 1:Q.n
            for i in 1:j
                if i == 1 && j == 2
                    x = el
                else
                    x = f(Q[i, j])
                end
                symmetric_setindex!(ret, x, i, j)
            end
        end
        return ret
    end
end

imag_map(n, i, j) = trimap(n, n) + trimap(i, j - 1)
imag_map(Q::VectorizedHermitianMatrix, i, j) = imag_map(Q.n, i, j)

function trimat(::Type{T}, f, n, σ) where {T}
    Q = Vector{T}(undef, N + trimap(n - 1, n - 1))
    for i in 1:n
        for j in 1:i
            x = f(σ[i], σ[j])
            Q[trimap(i, j)] = real(x)
            if i != j
                Q[imag_map(n, i, j)] = imag(x)
            end
        end
    end
    return VectorizedHermitianMatrix{T}(Q, n)
end

Base.size(Q::VectorizedHermitianMatrix) = (Q.n, Q.n)

"""
    symmetric_setindex!(Q::VectorizedHermitianMatrix, value, i::Integer, j::Integer)

Set `Q[i, j]` to the value `value` and `Q[j, i]` to the value `-value`.
"""
function symmetric_setindex!(Q::VectorizedHermitianMatrix, value, i::Integer, j::Integer)
    if i > j
        symmetric_setindex!(Q, conj(value), j, i)
    else
        Q.Q[trimap(i, j)] = real(value)
        if i == j
            if !iszero(imag(value))
                error("Cannot set diagonal entry ($i, $j) of hermitian matrix with a non-zero imaginary part $(imag(value))")
            end
        else
            Q.Q[imag_map(Q, i, j)] = imag(value)
        end
    end
end

function Base.getindex(Q::VectorizedHermitianMatrix, i::Integer, j::Integer)
    I, J = max(i, j), min(i, j)
    r = Q.Q[trimap(I, J)]
    if i == j
        return Complex(r)
    else
        c = Q.Q[imag_map(Q, I, J)]
        return Complex(r, i < j ? c : -c)
    end
end
Base.getindex(Q::VectorizedHermitianMatrix, I::Tuple) = Q[I...]
Base.getindex(Q::VectorizedHermitianMatrix, I::CartesianIndex) = Q[I.I]
