# Inspired from macaulaylab.net

struct MonomialDependence{B<:MB.AbstractPolynomialBasis}
    standard_monomials::B
    corners::B
    dependent_border::B
    independent_border::B
end

"""
    function standard_monomials_and_corners(
        null::MacaulayNullspace,
        rank_check,
    )

Computes the set of standard monomials using the *greedy sieve* algorithm
presented in [LLR08, Algorithm 1].

[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp.
*Semidefinite characterization and computation of zero-dimensional real radical ideals.*
Foundations of Computational Mathematics 8 (2008): 607-647.
"""
function standard_monomials_and_corners(
    null::MacaulayNullspace{T,MT,<:MB.MonomialBasis},
    rank_check,
) where {T,MT}
    monos = eltype(null.basis.monomials)[]
    corners = eltype(monos)[]
    rows = Int[]
    old_rank = 0
    for (k, mono) in enumerate(null.basis.monomials)
        if old_rank == size(null.matrix, 2)
            break
        end
        # This sieve of [LLR08, Algorithm 1] is a performance improvement but not only.
        # It also ensures that the standard monomials have the "staircase structure".
        if !any(Base.Fix2(MP.divides, mono), corners)
            new_rank =
                LinearAlgebra.rank(null.matrix[vcat(rows, k), :], rank_check)
            if new_rank < old_rank
                @warn(
                    "After adding rows, the rank dropped from `$old_rank` to `$new_rank`. Correcting the rank to `$old_rank` and continuing."
                )
                new_rank = old_rank
            elseif new_rank > old_rank
                push!(rows, k)
                push!(monos, mono)
            else
                push!(corner, mono)
            end
            old_rank = new_rank
        end
    end
    return monos, corners
end

function shift_nullspace(null::MacaulayNullspace, monos)
    S = null[monos]
    Sx = [null[monos.*shift] for shift in MP.variables(monos)]
    pS = LinearAlgebra.pinv(S.matrix)
    mult = SemialgebraicSets.MultiplicationMatrices([pS * S.matrix for S in Sx])
    return MultivariateMoments.SemialgebraicSets.solve(
        mult,
        MultivariateMoments.SemialgebraicSets.ReorderedSchurMultiplicationMatricesSolver{
            Float64,
        }(),
    )
end

function gap_zone_standard_monomials(monos, maxdegree)
    num = zeros(Int, maxdegree + 1)
    for mono in monos
        num[MP.maxdegree(mono)+1] += 1
    end
    i = findfirst(iszero, num)
    if isnothing(i)
        return
    end
    gap_size = something(
        findfirst(!iszero, @view(num[(i+1):end])),
        length(num) - i + 1,
    )
    num_affine = sum(view(num, 1:(i-1)))
    return num_affine, gap_size
end

"""
    struct ShiftNullspace{C<:RankCheck} <: MacaulayNullspaceSolver
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
struct ShiftNullspace{C<:RankCheck} <: MacaulayNullspaceSolver
    check::C
end
# Because the matrix is orthogonal, we know the SVD of the whole matrix is
# `ones(...)` so an `AbsoluteRankTol` would be fine here.
# However, since we also know that the first row (which correspond to the
# constant monomial) should be a standard monomial, `LeadingRelativeRankTol`
# ensures that we will take it.
ShiftNullspace() = ShiftNullspace(LeadingRelativeRankTol(1e-8))

function solve(null::MacaulayNullspace, shift::ShiftNullspace)
    d = MP.maxdegree(null.basis.monomials)
    monos, _ = standard_monomials_and_corners(null, shift.check)
    gap_zone = gap_zone_standard_monomials(monos, d)
    if isnothing(gap_zone)
        return
    end
    num_affine, gap_size = gap_zone
    if gap_size < 1
        return
    end
    if num_affine == length(monos)
        affine_monos = monos
        affine_null = null
    else
        affine_monos = monos[1:num_affine]
        @warn("Column compression not supported yet")
        return
    end

    # Solve the system:
    sols = shift_nullspace(affine_null, affine_monos)
    return ZeroDimensionalVariety(sols)
end
