"""
    quotient_basis(basis::MB.SubBasis{MB.Monomial}, domain)

Return a `MultivariateBases.QuotientBasis` representing the quotient
of the span of `basis` modulo the ideal of `domain`, which is either an
algebraic set, a polynomial ideal or a `SemialgebraicSets.TruncatedIdeal`. The
matrix is a sparse matrix with `float` of the coefficient type of `domain` as
element type.

If a Gröbner basis of the ideal is known, the matrix gives the remainders
of the elements of `basis` modulo the Gröbner basis. For a truncated ideal,
the polynomials whose coefficients are in the kernel of the matrix are the
polynomials in the row space of the Macaulay matrix of its degree. If `domain`
is an `AlgebraicSet` whose Gröbner basis is not known, it is not computed and
the ideal is instead truncated at the maximal degree of `basis`, see
`SemialgebraicSets.truncated`.
"""
function quotient_basis(
    basis::MB.SubBasis{MB.Monomial},
    domain::SS.AbstractAlgebraicSet,
)
    return quotient_basis(basis, SS.ideal(domain))
end

function quotient_basis(
    basis::MB.SubBasis{MB.Monomial},
    I::SS.AbstractPolynomialIdeal,
)
    T = float(MP.coefficient_type(typeof(I)))
    return MB.QuotientBasis(basis, _remainder_matrix(basis, I, T))
end

function quotient_basis(
    basis::MB.SubBasis{MB.Monomial},
    domain::SS.AlgebraicSet,
)
    I = SS.ideal(domain)
    if SS.is_gröbner_basis(I)
        return quotient_basis(basis, I)
    end
    return quotient_basis(basis, SS.truncated(I, _maxdegree(basis)))
end

_maxdegree(basis) = isempty(basis) ? 0 : MP.maxdegree(basis)

function quotient_basis(basis::MB.SubBasis{MB.Monomial}, I::SS.TruncatedIdeal)
    T = float(MP.coefficient_type(typeof(I)))
    W = _macaulay_zero_matrix(basis, I.generators, I.maxdegree, T)
    return MB.QuotientBasis(basis, W)
end

function MA.promote_operation(
    ::typeof(quotient_basis),
    ::Type{B},
    ::Type{D},
) where {T,I,B<:SA.AbstractBasis{T,I},D}
    M = SparseArrays.SparseMatrixCSC{float(MP.coefficient_type(D)),Int}
    return MB.QuotientBasis{T,I,B,M}
end

"""
    moment_vector(values::AbstractVector, basis::MB.QuotientBasis)

Return the moments of the elements of `basis.basis` given the values `values`
of the linear functional on the elements of `basis`. As the coefficients in
`basis` are `basis.matrix` times the coefficients in `basis.basis`, the
moments are obtained with the adjoint map `transpose(basis.matrix) * values`.
"""
function moment_vector(
    values::AbstractVector{T},
    basis::MB.QuotientBasis,
) where {T}
    _check_length(values, basis)
    return moment_vector(transpose(basis.matrix) * values, basis.basis)
end

# Sparse matrix `R` such that `R[i, j]` is the coefficient of the `i`th
# standard monomial in the remainder of the `j`th monomial of `basis` modulo
# `I`. Its number of nonzeros is the total number of terms of the remainders.
function _remainder_matrix(basis, I, ::Type{T}) where {T}
    if isempty(basis)
        return SparseArrays.spzeros(T, 0, 0)
    end
    rems = [rem(mono, I) for mono in MB.keys_as_monomials(basis)]
    standard = MP.monomial_vector(reduce(vcat, MP.monomials.(rems)))
    rows = Int[]
    cols = Int[]
    vals = T[]
    for (j, r) in enumerate(rems)
        for t in MP.terms(r)
            push!(rows, searchsortedfirst(standard, MP.monomial(t)))
            push!(cols, j)
            push!(vals, MP.coefficient(t))
        end
    end
    return SparseArrays.sparse(
        rows,
        cols,
        vals,
        length(standard),
        length(basis),
    )
end

# Sparse matrix `W` such that a polynomial `t` with support in `basis` is in
# the row space of the Macaulay matrix of degree `d` of `equalities` if and
# only if `W * coefficients(t, basis) = 0`.
function _macaulay_zero_matrix(basis, equalities, d, ::Type{T}) where {T}
    n = length(basis)
    if isempty(basis)
        return SparseArrays.spzeros(T, 0, 0)
    end
    monos = MB.keys_as_monomials(basis)
    # The shifts should be in all the variables and not only in the variables
    # of the equalities
    vars =
        MP.variables(MP.monomial_vector([monos; MP.monomials.(equalities)...]))
    M = MB.macaulay(equalities, d, vars)
    null = MacaulayNullspace(M)
    columns = MP.variables(null.basis)
    rows = Int[]
    found = Int[]
    others = Int[]
    for (i, mono) in enumerate(monos)
        j = SA.key_index(null.basis, MP.exponents(mono, columns))
        if isnothing(j)
            # The monomial is not in any row of the Macaulay matrix so it is
            # independent of the other ones
            push!(others, i)
        else
            push!(rows, j)
            push!(found, i)
        end
    end
    Z = null.matrix[rows, :]
    if isempty(Z)
        # Only monomials outside the Macaulay support survive in the quotient.
        return SparseArrays.sparse(
            eachindex(others),
            others,
            one(T),
            length(others),
            n,
        )
    end
    F = LinearAlgebra.svd(Z)
    rank_check = LeadingRelativeRankTol(Base.rtoldefault(T))
    r = rank_from_singular_values(F.S, rank_check)
    if r == length(found)
        # No monomial of `basis` is dependent so the polynomial must be zero
        return SparseArrays.sparse(one(T) * LinearAlgebra.I, n, n)
    end
    # The equations `transpose(F.U[:, 1:r]) * coefficients(t, basis) = 0` are
    # orthonormal but they are dense which is an issue for sparse solvers.
    # As any invertible recombination of the equations is equivalent, we
    # take the reduced row echelon form instead. The pivots are the smallest
    # independent monomials so if they coincide with the standard monomials
    # of a Gröbner basis, it coincides with the matrix of `_remainder_matrix`.
    E = Matrix(transpose(F.U[:, 1:r]))
    tol = Base.rtoldefault(T)
    RowEchelon.rref!(E, tol)
    # Remove the numerical noise so that the zeros are skipped
    E[abs.(E) .< tol] .= zero(T)
    # The pivot of each row of `W`, we sort the rows according to their
    # pivots so that the equations are in the order of `basis`
    pivots = [
        [found[findfirst(!iszero, view(E, k, :))] for k in 1:r]
        others
    ]
    W = zeros(T, length(pivots), n)
    W[1:r, found] = E
    for (k, i) in enumerate(others)
        W[r+k, i] = one(T)
    end
    return SparseArrays.sparse(W[sortperm(pivots), :])
end

"""
    BasisDependence{LinearDependence}(I::SemialgebraicSets.TruncatedIdeal, vars)

Return the [`BasisDependence`](@ref) of the monomials in the variables `vars`
that are columns of the Macaulay matrix of degree `I.maxdegree` of the
generators of `I`. A monomial `α` is dependent if there is a polynomial of
`I` whose leading monomial is `α`, i.e., `α` is congruent modulo `I` to a
combination of monomials smaller than `α`. Equivalently, the row of the null
space of the Macaulay matrix corresponding to `α` is a linear combination of
the rows of the monomials smaller than `α`. Monomials that are not columns
of the Macaulay matrix are independent.
"""
function BasisDependence{LinearDependence}(I::SS.TruncatedIdeal, vars)
    null = MacaulayNullspace(MB.macaulay(I.generators, I.maxdegree, vars))
    # We do not use the sieve of `BasisDependence{StaircaseDependence}` which
    # labels the multiples of dependent monomials as dependent. If `α` is
    # `d`-dependent, `x * α` may only be `(d + 1)`-dependent so this is not
    # valid for a truncated ideal, see
    # https://github.com/JuliaAlgebra/MultivariateMoments.jl/issues/106
    return BasisDependence{LinearDependence}(null)
end
