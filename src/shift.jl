# Inspired from macaulaylab.net

mutable struct RankDependence{T,MT<:AbstractMatrix{T},C}
    matrix::MT
    check::C
    independent_rows::Vector{Int}
    old_rank::Int
end

RankDependence(matrix, check) = RankDependence(matrix, check, Int[], 0)

function is_dependent!(r::RankDependence, row)
    if r.old_rank == size(r.matrix, 2)
        return true
    end
    rows = vcat(r.independent_rows, row)
    new_rank = LinearAlgebra.rank(r.matrix[rows, :], r.check)
    if new_rank < r.old_rank
        @warn(
            "After adding rows, the rank dropped from `$old_rank` to `$new_rank`. Correcting the rank to `$old_rank` and continuing."
        )
        new_rank = r.old_rank
    end
    independent = new_rank > r.old_rank
    r.old_rank = new_rank
    if independent
        r.independent_rows = rows
    end
    return !independent
end

function column_compression!(r::RankDependence, rows)
    if length(rows) < size(r.matrix, 2)
        # FIXME should we multiply by inv(Diagonal(s.S)) ?
        s = LinearAlgebra.svd(r.matrix[rows, :], full = true)
        r.matrix = r.matrix[dep_rows, :] * s.V
    end
end

function BasisDependence{LinearDependence}(null::MacaulayNullspace, rank_check::RankCheck)
    r = RankDependence(null.matrix, rank_check)
    return BasisDependence{LinearDependence}(r, null.basis)
end

function BasisDependence{StaircaseDependence}(null::MacaulayNullspace, rank_check::RankCheck)
    r = RankDependence(null.matrix, rank_check)
    return BasisDependence{StaircaseDependence}(r, null.basis)
end

function _indices(in::MB.MonomialBasis, from::MB.MonomialBasis)
    return Int[_index(in, mono) for mono in from.monomials]
end

function BorderBasis(d::BasisDependence{LinearDependence}, null::MacaulayNullspace)
    indep_rows = findall(d -> !is_dependent(d), d.dependence)
    dep_rows = findall(d -> is_dependent(d), d.dependence)
    @assert length(indep_rows) == size(null.matrix, 2)
    return BorderBasis(
        d,
        (null.matrix[dep_rows, :] / null.matrix[indep_rows, :])',
    )
end

function BorderBasis(d::StaircaseDependence, null::MacaulayNullspace)
    indep_rows = _indices(null.basis, standard_basis(d; trivial = false))
    dep_rows = _indices(null.basis, dependent_basis(d))
    U = null.matrix
    if length(indep_rows) < size(U, 2)
        U = LinearAlgebra.svd(U[indep_rows, :]).U
    end
    return BorderBasis(
        d,
        (null.matrix[dep_rows, :] / null.matrix[indep_rows, :])',
    )
end

function BorderBasis{D}(null::MacaulayNullspace, check::RankCheck) where {D}
    return BorderBasis(BasisDependence{D}(null, check), null)
end

"""
    struct ShiftNullspace{D,C<:RankCheck} <: MacaulayNullspaceSolver
        check::C
    end

From the [`MacaulayNullspace`](@ref), computes multiplication matrices
by exploiting the shift property of the rows [DBD12].
The rank check `check` is used to determine the standard monomials among
the row indices of the null space.

[DBD12] Dreesen, Philippe, Batselier, Kim, and De Moor, Bart.
*Back to the roots: Polynomial system solving, linear algebra, systems theory.*
IFAC Proceedings Volumes 45.16 (2012): 1203-1208.
"""
struct ShiftNullspace{D,C<:RankCheck} <: MacaulayNullspaceSolver
    check::C
end
# Because the matrix is orthogonal, we know the SVD of the whole matrix is
# `ones(...)` so an `AbsoluteRankTol` would be fine here.
# However, since we also know that the first row (which correspond to the
# constant monomial) should be a standard monomial, `LeadingRelativeRankTol`
# ensures that we will take it.
function ShiftNullspace{D}(check::RankCheck) where {D}
    return ShiftNullspace{D,typeof(check)}(check)
end
ShiftNullspace{D}() where {D} = ShiftNullspace{D}(LeadingRelativeRankTol(1e-8))
ShiftNullspace(args...) = ShiftNullspace{StaircaseDependence}(args...)

function border_basis_and_solver(
    null::MacaulayNullspace,
    shift::ShiftNullspace{D},
) where {D}
    return BorderBasis{D}(null, shift.check), nothing
end
