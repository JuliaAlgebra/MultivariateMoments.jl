using SemialgebraicSets

abstract type AbstractMomentMatrix{T,B<:SA.ExplicitBasis} <:
              AbstractMeasureLike{T} end

"""
    mutable struct MomentMatrix{T,B<:SA.ExplicitBasis,MT<:AbstractMatrix{T}} <: AbstractMeasureLike{T}
        Q::MT
        basis::B
        support::Union{Nothing, AlgebraicSet}
    end

Measure ``\\nu`` represented by the moments of the monomial matrix ``x x^\\top`` in the symmetric matrix `Q`.
The set of points that are zeros of all the polynomials ``p`` such that ``\\mathbb{E}_{\\nu}[p] = 0`` is stored in the field `support` when it is computed.
"""
mutable struct MomentMatrix{T,B<:SA.ExplicitBasis,MT<:AbstractMatrix{T}} <:
               AbstractMomentMatrix{T,B}
    Q::MT
    basis::B
    support::Union{Nothing,AbstractAlgebraicSet}
end
function MomentMatrix{T,B,MT}(Q::MT, basis::SA.ExplicitBasis) where {T,B,MT}
    return MomentMatrix{T,B,MT}(Q, basis, nothing)
end
function MomentMatrix{T,B}(
    Q::AbstractMatrix{T},
    basis::SA.ExplicitBasis,
) where {T,B}
    return MomentMatrix{T,B,typeof(Q)}(Q, basis)
end
function MomentMatrix(Q::SymMatrix{T}, basis::SA.ExplicitBasis) where {T}
    return MomentMatrix{T,typeof(basis)}(Q, basis)
end

MP.variables(μ::MomentMatrix) = MP.variables(μ.basis)
MP.nvariables(μ::MomentMatrix) = MP.nvariables(μ.basis)

function MomentMatrix{T}(
    f::Function,
    basis::SA.ExplicitBasis,
    σ = 1:length(basis),
) where {T}
    return MomentMatrix(
        vectorized_symmetric_matrix(T, f, length(basis), σ),
        basis,
    )
end
function MomentMatrix{T}(f::Function, monos::AbstractVector) where {T}
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return MomentMatrix{T}(f, MB.SubBasis{MB.Monomial}(sorted_monos), σ)
end

function show_basis_indexed_matrix(io::IO, A, pre = "")
    println(io, " with row/column basis:")
    println(io, pre, " ", A.basis)
    print(io, pre, "And entries in a ", summary(A.Q))
    isempty(A.Q) && return
    println(io, ":")
    Base.print_matrix(io, A.Q, pre * " ")
    return
end

function Base.show(io::IO, M::MomentMatrix)
    print(io, "MomentMatrix")
    return show_basis_indexed_matrix(io, M)
end

"""
    moment_matrix(μ::MomentVector, x)

Creates a matrix the moment matrix for the moment matrix  ``x x^\\top`` using the moments of `μ`.
"""
function moment_matrix(μ::MomentVector{T}, basis) where {T}
    return MomentMatrix{T}(
        (i, j) -> moment_value(μ, basis[i] * basis[j]),
        basis,
    )
end

function MomentMatrix(
    Q::AbstractMatrix{T},
    basis::SA.ExplicitBasis,
    σ,
) where {T}
    return MomentMatrix{T}((i, j) -> Q[σ[i], σ[j]], basis)
end
function MomentMatrix(Q::AbstractMatrix, monos::AbstractVector)
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return MomentMatrix(Q, MB.SubBasis{MB.Monomial}(sorted_monos), σ)
end
moment_matrix(Q::AbstractMatrix, monos) = MomentMatrix(Q, monos)

value_matrix(μ::MomentMatrix) = Matrix(μ.Q)

function vectorized_basis(
    ν::MomentMatrix{T,<:MB.SubBasis{MB.Monomial}},
) where {T}
    monos = ν.basis.monomials
    n = length(monos)
    # We don't wrap in `MB.SubBasis` as we don't want the monomials
    # to be `sort`ed and `uniq`ed.
    return [monos[i] * monos[j] for j in 1:n for i in 1:j]
end

function moment_vector(ν::MomentMatrix; kws...)
    n = length(ν.basis)
    return moment_vector(ν.Q.Q, vectorized_basis(ν); kws...)
end

struct BlockDiagonalMomentMatrix{T,B<:SA.ExplicitBasis,MT} <:
       AbstractMomentMatrix{T,B}
    blocks::Vector{MomentMatrix{T,B,MT}}
end

function block_diagonal(blocks::Vector{<:MomentMatrix})
    return BlockDiagonalMomentMatrix(blocks)
end

function show_basis_indexed_blocks(io::IO, blocks)
    print(io, " with $(length(blocks)) blocks:")
    nd = ndigits(length(blocks))
    for i in eachindex(blocks)
        println(io)
        print(io, "[", " "^(nd - ndigits(i)), i, "] Block")
        show_basis_indexed_matrix(io, blocks[i], " "^(nd + 3))
    end
    return
end

function Base.show(io::IO, M::BlockDiagonalMomentMatrix)
    print(io, "BlockDiagonalMomentMatrix")
    show_basis_indexed_blocks(io, M.blocks)
    return
end
