function build_system(U::AbstractMatrix, basis::MB.MonomialBasis, ztol, args...)
    # System is
    # y = [U 0] * y
    # where y = x[end:-1:1]
    # which is
    # y = U * β
    m = length(basis)
    r = size(U, 1)
    pivots = [findfirst(j -> U[i, j] != 0, 1:m) for i in 1:r]
    if any(isnothing, pivots)
        keep = map(!isnothing, pivots)
        pivots = pivots[keep]
        r = length(pivots)
        U = U[keep, :]
    end
    monos = basis.monomials
    β = monos[pivots]
    system = [MA.operate(dot, β, U[:, i]) - monos[i] for i in eachindex(monos)]
    filter!(!MP.isconstant, system)
    # Type instability here :(
    if MP.mindegree(monos) == MP.maxdegree(monos) # Homogeneous
        projective_algebraic_set(system, Buchberger(ztol), args...)
    else
        algebraic_set(system, Buchberger(ztol), args...)
    end
end

"""
    struct Echelon
    end

Given a [`MacaulayNullspace`](@ref), computes its echelon form
(corresponding to the *Canonical Null Space* of [D13]) with
Gaussian elimination.
From this echelon form, the left null space can easily be computed
using using [HL05, (8)].
This left null space forms a system of polynomial equations.

!!! note
    In the context of [`compute_support!`](@ref), if the
    [`MacaulayNullspace`](@ref) was obtained using [`SVDLDLT`](@ref), the
    left null space can easily be obtained from the singular vectors
    corresponding to the negligeable singular values that were removed.
    However, as mentioned [LLR08, Section 4.4.5], these will usually give an
    overdetermined bases.
    As shown in [S04, Section 10.8.2], it is desirable to avoid overdetermined
    bases because it could lead to inconsistencies in the basis for numerical
    reasons.
    For this reason, this method computes the left null space from the
    echelon form of the significant singular values instead.

Let `B` be the set of monomials corresponding to the rows of the pivots of this
echelon form.  If the moment matrix satisfies the *flat extension* property
described in [L09, Section 5.3] then all monomials of the *border* of `B` (as
defined in [LLR08, (2.3)]) will correspond to a row of of the matrix.
In that case, the polynomial of the system obtained by [HL05, (8)] form a
*rewriting family* for `B` [L09, (2.16)] a.k.a. a `B`*-border prebasis* (as
defined in [LLR08, (2.4)]).
Therefore, they can be used to compute multiplication matrices.

[HL05] Henrion, D. & Lasserre, J-B.
*Detecting Global Optimality and Extracting Solutions of GloptiPoly*
2005

[D13] Dreesen, Philippe.
*Back to the Roots: Polynomial System Solving Using Linear Algebra*
Ph.D. thesis (2013)

[L09] Laurent, Monique.
*Sums of squares, moment matrices and optimization over polynomials.*
Emerging applications of algebraic geometry (2009): 157-270.

[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp.
*Semidefinite characterization and computation of zero-dimensional real radical ideals.*
Foundations of Computational Mathematics 8 (2008): 607-647.

[S04] Stetter, Hans J.
*Numerical polynomial algebra.*
Society for Industrial and Applied Mathematics, 2004.
"""
struct Echelon end

import RowEchelon

function solve(null::MacaulayNullspace, ::Echelon, args...)
    # If M is multiplied by λ, W is multiplied by √λ
    # so we take √||M|| = √nM
    Z = Matrix(null.matrix')
    RowEchelon.rref!(Z, null.accuracy / sqrt(size(Z, 2)))
    #r, vals = solve_system(U', μ.x)
    # TODO determine what is better between rank_check and sqrt(rank_check) here
    return build_system(Z, null.basis, √null.accuracy, args...)
end
