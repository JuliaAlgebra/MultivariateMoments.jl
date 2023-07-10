# Inspired from macaulaylab.net

# Cannot call it as exported symbol `standard_monomials` as it would
# collide with `SemialgebraicSets.standard_monomials` and cannot add a method
# as it would be type piracy
# TODO implement sieve
function _standard_monomials(Z, tol = 1e-10)
    list = Int[]
    old_rank = 0
    for k in axes(Z, 1)
        new_rank = LinearAlgebra.rank(Z[1:k, :], tol)
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

function shiftnullspace(Z, dgap, srows, monos)
    S = Z[srows, :]
    Sx = [
        Z[[findfirst(isequal(monos[row] * x), monos) for row in srows], :]
        for x in MP.variables(monos)
    ]
    pS = LinearAlgebra.pinv(S)
    mult = SemialgebraicSets.MultiplicationMatrices([pS * S for S in Sx])
    return MultivariateMoments.SemialgebraicSets.solve(
        mult,
        MultivariateMoments.SemialgebraicSets.ReorderedSchurMultiplicationMatricesSolver{
            Float64,
        }(),
    )

    return MultivariateMoments.solve
    eig = [LinearAlgebra.eigen(S).values for S in Sx]
    return [[eig[i][j] for i in eachindex(eig)] for j in eachindex(eig[1])]
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
struct ShiftNullspace end

function solve(null::MacaulayNullspace, ::ShiftNullspace)
    monos = null.basis.monomials
    Z = null.matrix
    d = MP.maxdegree(monos)
    srows = _standard_monomials(Z)
    gap_zone = gap_zone_standard_monomials(monos[srows], d)
    if isnothing(gap_zone)
        return
    end
    dgap, ma, gapsize = gap_zone
    srows = srows[1:ma]
    if gapsize < 1
        return
    end
    mb = size(Z, 2)
    if mb == ma
        W = Z
    else
        @warn("Column compression not supported yet")
        return
    end

    # Solve the system:
    sols = shiftnullspace(W, dgap, srows, monos)
    return ZeroDimensionalVariety(sols)
end
