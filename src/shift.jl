# Inspired from macaulaylab.net

# Cannot call it as exported symbol `standard_monomials` as it would
# collide with `SemialgebraicSets.standard_monomials` and cannot add a method
# as it would be type piracy
# TODO implement sieve
function _standard_monomials(Z, rank_check)
    list = Int[]
    old_rank = 0
    for k in axes(Z, 1)
        new_rank = LinearAlgebra.rank(Z[1:k, :], rank_check)
        if new_rank > old_rank
            push!(list, k)
        end
        old_rank = new_rank
        if new_rank == size(Z, 2)
            break
        end
    end
    return list
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
    dgap = i - 2
    gapsize = something(
        findfirst(!iszero, @view(num[(i+1):end])),
        length(num) - i + 1,
    )
    ma = sum(view(num, 1:(i-1)))
    return dgap, ma, gapsize
end

"""
    struct ShiftNullspace end

From the [`MacaulayNullspace`](@ref), computes multiplication matrices
by exploiting the shift property of the rows [DBD12].

[DBD12] Dreesen, Philippe, Batselier, Kim, and De Moor, Bart.
*Back to the roots: Polynomial system solving, linear algebra, systems theory.*
IFAC Proceedings Volumes 45.16 (2012): 1203-1208.
"""
struct ShiftNullspace{C<:RankCheck}
    check::C
end
# Because the matrix is orthogonal, we know the SVD of the whole matrix is
# `ones(...)` so an `AbsoluteRankTol` would be fine here.
# However, since we also know that the first row (which correspond to the
# constant monomial) should be a standard monomial, `LeadingRelativeRankTol`
# ensures that we will take it.
ShiftNullspace() = ShiftNullspace(LeadingRelativeRankTol(1e-8))

function solve(null::MacaulayNullspace, shift::ShiftNullspace)
    Z = null.matrix
    d = MP.maxdegree(null.basis.monomials)
    srows = _standard_monomials(Z, shift.check)
    monos = null.basis.monomials[srows]
    gap_zone = gap_zone_standard_monomials(monos, d)
    if isnothing(gap_zone)
        return
    end
    dgap, ma, gapsize = gap_zone
    if gapsize < 1
        return
    end
    mb = size(Z, 2)
    if mb == ma
        affine_monos = monos
        affine_null = null
    else
        affine_monos = monos[1:ma]
        @warn("Column compression not supported yet")
        return
    end

    # Solve the system:
    sols = shift_nullspace(affine_null, affine_monos)
    return ZeroDimensionalVariety(sols)
end
