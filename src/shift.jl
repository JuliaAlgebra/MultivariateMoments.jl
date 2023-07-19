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
        return false
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
    return independent
end

function AnyDependence(null::MacaulayNullspace, rank_check::RankCheck)
    r = RankDependence(null.matrix, rank_check)
    return AnyDependence(Base.Fix1(is_dependent!, r), null.basis)
end

function StaircaseDependence(null::MacaulayNullspace, rank_check::RankCheck)
    r = RankDependence(null.matrix, rank_check)
    return StaircaseDependence(Base.Fix1(is_dependent!, r), null.basis)
end

"""
"""
function column_compression end

function _indices_or_ignore(in::MB.MonomialBasis, from::MB.MonomialBasis)
    indices = Int[]
    for mono in from.monomials
        i = _index(in, mono)
        if !isnothing(i)
            push!(indices, i)
        end
    end
    return indices
end

function _indices(in::MB.MonomialBasis, from::MB.MonomialBasis)
    return Int[_index(in, mono) for mono in from.monomials]
end

function BorderBasis(d::AnyDependence, null::MacaulayNullspace)
    indep_rows = _indices_or_ignore(null.basis, d.independent)
    dep_rows = _indices(null.basis, d.dependent)
    @assert length(indep_rows) == size(null.matrix, 2)
    return BorderBasis(d, null.matrix[dep_rows, :] / null.matrix[indep_rows, :])
end

function BorderBasis(d::StaircaseDependence, null::MacaulayNullspace)
    indep_rows = _indices_or_ignore(null.basis, d.independent)
    dep = MB.merge_bases(d.corners, d.dependent_border)
    dep_rows = _indices(null.basis, dep)
    if length(indep_rows) < size(null.matrix)
        error("Column compression not supported yet")
    end
    return BorderBasis(d, null.matrix[dep_rows, :] / null.matrix[indep_rows, :])
end

function BorderBasis{D}(null::MacaulayNullspace, check::RankCheck) where {D}
    return BorderBasis(D(null, check), null)
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
