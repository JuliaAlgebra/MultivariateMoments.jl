"""
    struct MacaulayNullspace{T,MT<:AbstractMatrix{T},BT}
        matrix::MT
        basis::BT
        accuracy::T
    end

This matrix with rows indexed by `basis` can either be interpreted as
the right null space of a Macaulay matrix with columns indexed by `basis` (resp. or
the image space of a moment matrix with rows and columns indexed by `basis`).
The value of `matrix[i, j]` should be interpreted as the value of the `i`th
element of `basis` for the `j`th generator of the null space (resp. image) space.
"""
struct MacaulayNullspace{T,MT<:AbstractMatrix{T},BT<:SA.ExplicitBasis}
    matrix::MT
    basis::BT
    accuracy::T
end
function MacaulayNullspace(matrix::AbstractMatrix{T}, basis) where {T}
    return MacaulayNullspace(matrix, basis, Base.rtoldefault(T))
end

function Base.getindex(
    null::MacaulayNullspace{T,MT,<:MB.SubBasis{B}},
    monos,
) where {T,MT,B}
    I = _index.(Ref(null.basis), monos)
    return MacaulayNullspace(
        null.matrix[I, :],
        MB.SubBasis{B}(null.basis.monomials[I]),
        null.accuracy,
    )
end

function MacaulayNullspace(
    ν::MomentMatrix,
    rank_check::RankCheck,
    ldlt::LowRankLDLTAlgorithm = SVDLDLT(),
)
    M = value_matrix(ν)
    chol = low_rank_ldlt(M, ldlt, rank_check)
    @assert size(chol.L, 1) == LinearAlgebra.checksquare(M)
    return MacaulayNullspace(chol.L, ν.basis, accuracy(chol))
end

abstract type MacaulayNullspaceSolver end

function solve(null::MacaulayNullspace, solver::MacaulayNullspaceSolver)
    border, solver = border_basis_and_solver(null, solver)
    return solve(border, _some_args(solver)...)
end

"""
    struct ImageSpaceSolver{A<:LowRankLDLTAlgorithm,S<:MacaulayNullspaceSolver}
        ldlt::A
        null::S
    end

Computes the image space of the moment matrix using `ldlt` and then solve it
with `null`.
"""
struct ImageSpaceSolver{A<:LowRankLDLTAlgorithm,S<:MacaulayNullspaceSolver}
    ldlt::A
    null::S
end

function compute_support!(
    ν::MomentMatrix,
    rank_check::RankCheck,
    solver::ImageSpaceSolver,
)
    ν.support = solve(MacaulayNullspace(ν, rank_check, solver.ldlt), solver.null)
    return
end

function compute_support!(
    ν::MomentMatrix,
    rank_check::RankCheck,
    ldlt::LowRankLDLTAlgorithm,
    null::MacaulayNullspaceSolver,
)
    return compute_support!(ν, rank_check, ImageSpaceSolver(ldlt, null))
end

function compute_support!(
    ν::MomentMatrix,
    rank_check::RankCheck,
    ldlt::LowRankLDLTAlgorithm,
    args...,
)
    return compute_support!(ν, rank_check, ldlt, Echelon(args...))
end

function compute_support!(μ::MomentMatrix, rank_check::RankCheck, args...)
    return compute_support!(
        μ::MomentMatrix,
        rank_check::RankCheck,
        SVDLDLT(),
        args...,
    )
end
